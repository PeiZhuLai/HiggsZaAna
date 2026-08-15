"""PyTorch reconstruction of the MLPhoton classifier and mass regressor.

Layer-for-layer identical to the deployed ONNX graphs (pytorch 1.7, opset 9), so
that:
  * the existing weights can be loaded back in and reproduce the current models
    bit for bit -- verified by verify_models.py, the model-level analogue of the
    preprocessing closure test;
  * a retrained model exports to a graph CMSSW's MLPhotonProducer can load
    without any C++ change.

Shapes, read out of the ONNX attributes rather than guessed:

  regressor   30 -conv7p2-> 28 -conv5p2-> 28 -pool2-> 14 -conv4p2-> 15
              -conv3p1-> 15 -pool2-> 7   =>  16*7*7 = 784
  classifier  30 -conv7p2-> 28 -conv5p2-> 28 -pool2-> 14 -conv4p2-> 15
              -conv3p2-> 17 -pool2-> 8   =>  64*8*8 = 4096

Note conv4 is padded by 1 in the regressor but by 2 in the classifier. That
asymmetry is in the deployed graphs; it is reproduced deliberately.

Two more inherited details that must not be "tidied up":
  * there is NO activation between conv4 and the second max-pool;
  * the regressor ends in a LeakyReLU on the scalar output, so it can emit small
    negative m/E values. The deployed model does this, and the analysis has
    always seen it.
"""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

ISIZE = 30
LEAKY_SLOPE = 0.01


class MLPhotonClassifier(nn.Module):
    """3-class cluster classifier: index 0 = mono, 1 = di, 2 = hadronic.

    Input is the RAW (unnormalised) energy image -- see cluster_py.
    Output is log-softmax, matching the deployed graph; MLPhotonProducer then
    exponentiates and renormalises, which is a no-op on top of log-softmax.
    """

    def __init__(self, nf: int = 64, n_classes: int = 3):
        super().__init__()
        self.conv1 = nn.Conv2d(1, nf, 7, padding=2)
        self.conv2 = nn.Conv2d(nf, nf, 5, padding=2)
        self.conv3 = nn.Conv2d(nf, nf, 4, padding=2)
        self.conv4 = nn.Conv2d(nf, nf, 3, padding=2)
        self.fc1 = nn.Linear(nf * 8 * 8, n_classes)

    def forward(self, img):
        x = F.relu(self.conv1(img))
        x = F.relu(self.conv2(x))
        x = F.max_pool2d(x, 2, 2)
        x = F.relu(self.conv3(x))
        x = self.conv4(x)          # no activation here, as in the ONNX graph
        x = F.max_pool2d(x, 2, 2)
        x = x.reshape(x.size(0), -1)
        x = self.fc1(x)
        return F.log_softmax(x, dim=1)


class MLPhotonRegressor(nn.Module):
    """Regresses m/E from the SUM-NORMALISED image plus the cluster eta."""

    def __init__(self, nf: int = 16):
        super().__init__()
        self.conv1 = nn.Conv2d(1, nf, 7, padding=2)
        self.conv2 = nn.Conv2d(nf, nf, 5, padding=2)
        self.conv3 = nn.Conv2d(nf, nf, 4, padding=2)
        self.conv4 = nn.Conv2d(nf, nf, 3, padding=1)   # pad 1 here, 2 in the classifier
        self.fc1 = nn.Linear(nf * 7 * 7, 64)
        self.fc2 = nn.Linear(64 + 1, 16)
        self.fc3 = nn.Linear(16, 1)

    def forward(self, img, eta):
        x = F.relu(self.conv1(img))
        x = F.relu(self.conv2(x))
        x = F.max_pool2d(x, 2, 2)
        x = F.relu(self.conv3(x))
        x = self.conv4(x)          # no activation here
        x = F.max_pool2d(x, 2, 2)
        x = x.reshape(x.size(0), -1)
        x = F.leaky_relu(self.fc1(x), LEAKY_SLOPE)
        x = torch.cat([x, eta.reshape(x.size(0), 1)], dim=1)
        x = F.leaky_relu(self.fc2(x), LEAKY_SLOPE)
        x = self.fc3(x).squeeze(1)
        return F.leaky_relu(x, LEAKY_SLOPE)


class MLPhotonRegressorLog(MLPhotonRegressor):
    """Same trunk, but predicts m/E through a standardised LOG target.

    Why not keep the deployed head: the target spans 1e-4 to 0.9. Training the
    deployed head on a relative loss collapses -- with
    L = ((pred-target)/target)^2, predicting zero costs exactly 1.0 per sample
    while any overshoot costs much more, so the optimiser parks the output at
    zero. Measured: from scratch it reaches bias = -0.001 and val loss 1.00,
    which is WORSE than the untrained deployed model (0.63). That is a property
    of the loss, not of the data.

    So the network emits a standardised log instead,

        raw   = trunk(image, eta)                       (unbounded, linear)
        m/E   = exp(raw * sigma + mu)

    and is trained with plain MSE on `raw` against (log(target) - mu)/sigma.
    The output is positive by construction, every mass point contributes on the
    same scale, and there is no zero-collapse basin.

    The exported graph still ends in a single scalar m/E, which is all
    MLPhotonProducer consumes -- it reads `regress_outputs.at(0)` and multiplies
    by the cluster energy. The extra Mul/Add/Exp nodes need no C++ change.

    mu/sigma are buffers, so they travel with the state_dict and get folded into
    the ONNX graph on export. Set them from the training sample before fitting.
    """

    def __init__(self, nf: int = 16, mu: float = 0.0, sigma: float = 1.0):
        super().__init__(nf=nf)
        self.register_buffer("mu", torch.tensor(float(mu)))
        self.register_buffer("sigma", torch.tensor(float(sigma)))

    def trunk(self, img, eta):
        """Standardised log(m/E). This is what the loss is computed on."""
        x = F.relu(self.conv1(img))
        x = F.relu(self.conv2(x))
        x = F.max_pool2d(x, 2, 2)
        x = F.relu(self.conv3(x))
        x = self.conv4(x)
        x = F.max_pool2d(x, 2, 2)
        x = x.reshape(x.size(0), -1)
        x = F.leaky_relu(self.fc1(x), LEAKY_SLOPE)
        x = torch.cat([x, eta.reshape(x.size(0), 1)], dim=1)
        x = F.leaky_relu(self.fc2(x), LEAKY_SLOPE)
        return self.fc3(x).squeeze(1)      # linear: no activation on the head

    def forward(self, img, eta):
        """m/E, positive by construction -- this is what gets exported."""
        return torch.exp(self.trunk(img, eta) * self.sigma + self.mu)

    def set_target_scale(self, log_targets) -> None:
        with torch.no_grad():
            self.mu.fill_(float(np.mean(log_targets)))
            self.sigma.fill_(float(np.std(log_targets)))


# --------------------------------------------------------------------------
# ONNX weight transfer
# --------------------------------------------------------------------------
def load_onnx_weights(model: nn.Module, onnx_path: str, strict: bool = True) -> list:
    """Copy initializers out of an ONNX file into a matching torch module.

    The ONNX initializer names are already 'conv1.weight', 'fc2.bias', ... so
    they map straight onto the torch parameter names. Returns the list of names
    that were transferred.
    """
    import onnx
    from onnx import numpy_helper

    graph = onnx.load(onnx_path).graph
    weights = {init.name: numpy_helper.to_array(init) for init in graph.initializer}

    own = dict(model.named_parameters())
    missing = [k for k in own if k not in weights]
    unused = [k for k in weights if k not in own]
    if strict and (missing or unused):
        raise RuntimeError(
            f"weight name mismatch for {onnx_path}\n  missing in onnx: {missing}"
            f"\n  unused from onnx: {unused}")

    transferred = []
    with torch.no_grad():
        for name, param in own.items():
            if name not in weights:
                continue
            arr = weights[name]
            if tuple(arr.shape) != tuple(param.shape):
                raise RuntimeError(
                    f"shape mismatch for {name}: onnx {arr.shape} vs torch "
                    f"{tuple(param.shape)}")
            param.copy_(torch.from_numpy(np.array(arr)))
            transferred.append(name)
    return transferred


def build_from_onnx(kind: str, onnx_path: str) -> nn.Module:
    """kind: 'classifier' or 'regressor'."""
    model = MLPhotonClassifier() if kind == "classifier" else MLPhotonRegressor()
    load_onnx_weights(model, onnx_path)
    model.eval()
    return model


# --------------------------------------------------------------------------
# ONNX export
# --------------------------------------------------------------------------
def _export(model, args_tuple, path, input_names, opset):
    """Always take the TorchScript exporter.

    torch >= 2.9's `torch.onnx.export` defaults to the dynamo path, which pulls
    in `onnxscript` at import time -- not installed in `hzgdna`, and we do not
    add packages to shared envs. `torch.onnx.utils.export` is the classic path
    and needs nothing extra. It also emits exactly the opset-9 graph shape the
    deployed models use.

    The "squeeze on dimension 1 / use opset 11 for dynamic shapes" warning is
    expected and harmless: MLPhotonProducer runs one cluster at a time, so the
    batch dimension is always 1.
    """
    model.eval()
    torch.onnx.utils.export(model, args_tuple, path, opset_version=opset,
                            input_names=input_names, output_names=["output"],
                            do_constant_folding=True)


def export_classifier(model: nn.Module, path: str, opset: int = 9) -> None:
    """Export with the input/output names MLPhotonProducer looks up ('img')."""
    _export(model, (torch.zeros(1, 1, ISIZE, ISIZE),), path, ["img"], opset)


def export_regressor(model: nn.Module, path: str, opset: int = 9) -> None:
    """Export with names 'img' and 'eta' -- MLPhotonProducer passes them in that
    order and looks them up by name, so both must match exactly."""
    _export(model, (torch.zeros(1, 1, ISIZE, ISIZE), torch.zeros(1, 1)),
            path, ["img", "eta"], opset)

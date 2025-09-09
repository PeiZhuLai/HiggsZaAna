#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
并行同步 Parquet 与 cutflow .out（Python 版，抗抖增强，兼容 Py3.8/3.9）
用法：
  python 2_move_files.py
  python 2_move_files.py --jobs 3 --logdir ./logs_sync
"""

import argparse
import asyncio
import os
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Tuple, Callable

# ===== 默认路径（可用参数覆盖）=====
SRC_ROOT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet"
DST_PARQ = "/eos/home-p/pelai/HZa/Parquet/NanoV12/run3"
DST_OUT  = "/afs/cern.ch/work/p/pelai/HZa/Output/cutflow_outfile/run3"

# rsync 旗标：更稳健的组合（遇到抖动可续传+校验；设置超时避免永久挂住）
RSYNC_FLAGS: List[str] = ["-a", "--info=progress2", "--partial", "--append-verify", "--timeout=600"]

# SUBDIRS: List[str] = ["Bkg_MC", "Data", "Sig_MC"]
SUBDIRS: List[str] = ["Bkg_MC"]

# 用 stdbuf 取消缓冲（几乎所有 Linux 都有）
STDBUF_PREFIX: List[str] = ["stdbuf", "-oL", "-eL"]


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


async def _stream_to_log(stream: asyncio.StreamReader, logf):
    """
    按塊讀取 stdout/stderr，避免 rsync 輸出超長單行導致 asyncio.readline() 的上限報錯。
    同時把進度輸出的 '\r' 轉為 '\n'，讓日誌可閱讀。
    """
    CHUNK = 16384  # 16KB 一塊，按需可調
    while True:
        chunk = await stream.read(CHUNK)
        if not chunk:
            break
        text = chunk.decode(errors="replace")
        # rsync 的進度列常用 '\r' 刷新，寫到檔案時轉成換行更清楚
        text = text.replace("\r", "\n")
        logf.write(text)
        logf.flush()
        # 控制台打一顆點，可自行移除
        try:
            sys.stdout.write(".")
            sys.stdout.flush()
        except Exception:
            pass


async def run_cmd_with_live_log(cmd: List[str], log_path: Path, stdin_bytes: Optional[bytes] = None) -> int:
    """
    启动子进程，实时把 stdout/stderr 写入 log_path。
    若提供 stdin_bytes，则写入进程 stdin（用于 rsync --files-from=-）。
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        stdin=asyncio.subprocess.PIPE if stdin_bytes is not None else None,
        start_new_session=True,   # 与控制终端会话分离，避免 SIGHUP 直接杀进程
    )

    with open(log_path, "a", buffering=1) as logf:
        logf.write("[{}] START: {}\n".format(ts(), " ".join(cmd)))
        logf.flush()

        async def _feed_stdin():
            if stdin_bytes is None or proc.stdin is None:
                return
            try:
                proc.stdin.write(stdin_bytes)
                await proc.stdin.drain()
            finally:
                try:
                    proc.stdin.close()
                except Exception:
                    pass

        feeder = asyncio.create_task(_feed_stdin())

        tasks: List[asyncio.Task] = []
        if proc.stdout:
            tasks.append(asyncio.create_task(_stream_to_log(proc.stdout, logf)))
        if proc.stderr:
            tasks.append(asyncio.create_task(_stream_to_log(proc.stderr, logf)))

        await asyncio.gather(*tasks)
        await feeder
        rc = await proc.wait()
        logf.write("\n[{}] DONE (rc={}): {}\n".format(ts(), rc, " ".join(cmd)))
        logf.flush()

    return rc


async def run_with_retry(cmd_builder: Callable[[], List[str]],
                         log_path: Path,
                         stdin_bytes: Optional[bytes] = None,
                         retries: int = 5,
                         retry_on: Tuple[int, ...] = (19, 20, 30)) -> int:
    """
    对 rsync 常见的信号类退出码做重试：
      19: SIGUSR1（伴随 20）
      20: SIGINT/SIGTERM/SIGHUP（常见于会话/IO 抖动）
      30: 超时
    指数退避：5s -> 10s -> 20s -> 40s -> 60s（封顶）
    """
    backoff = 5
    last_rc = 0
    for i in range(retries):
        cmd = cmd_builder()
        rc = await run_cmd_with_live_log(cmd, log_path, stdin_bytes=stdin_bytes)
        last_rc = rc
        if rc not in retry_on:
            return rc
        with open(log_path, "a", buffering=1) as f:
            f.write("[{}] WARN: rc={} matched {}, retry {}/{} after {}s\n"
                    .format(ts(), rc, retry_on, i+1, retries, backoff))
        await asyncio.sleep(backoff)
        backoff = min(backoff * 2, 60)
    return last_rc


def build_out_files_stdin_list(src_dir: Path) -> bytes:
    """
    收集 src_dir 下所有以 .out 结尾的相对路径，用 NUL 分隔，返回 bytes。
    等价于：find SRC -type f -name "*.out" -printf "%P\\0"
    """
    items: List[bytes] = []
    base = src_dir.resolve()
    for root, _, files in os.walk(base):
        for name in files:
            if name.endswith(".out"):
                full = Path(root) / name
                rel = full.relative_to(base)
                items.append((str(rel)).encode("utf-8") + b"\x00")
    return b"".join(items)


async def run_parquet_sync(src_root: str, dst_parq: str, subdir: str, logdir: Path) -> int:
    """Parquet 全量同步"""
    src = "{}/{}/".format(src_root.rstrip("/"), subdir)
    dst = "{}/{}/".format(dst_parq.rstrip("/"), subdir)
    Path(dst).mkdir(parents=True, exist_ok=True)
    log_path = logdir / ("parquet_{}.log".format(subdir.replace("/", "_")))

    def cmd_builder() -> List[str]:
        return STDBUF_PREFIX + ["rsync"] + RSYNC_FLAGS + [src, dst]

    return await run_with_retry(cmd_builder, log_path)


async def run_out_sync(src_root: str, dst_out: str, subdir: str, logdir: Path) -> int:
    """仅同步 .out 文件"""
    src_dir = Path("{}/{}".format(src_root.rstrip("/"), subdir))
    dst_dir = Path("{}/{}".format(dst_out.rstrip("/"), subdir))
    dst_dir.mkdir(parents=True, exist_ok=True)

    log_path = logdir / ("out_{}.log".format(subdir.replace("/", "_")))
    files_bytes = build_out_files_stdin_list(src_dir)
    if not files_bytes:
        with open(log_path, "a", buffering=1) as f:
            f.write("[{}] INFO: no .out under {}\n".format(ts(), src_dir))
        return 0

    def cmd_builder() -> List[str]:
        return STDBUF_PREFIX + ["rsync"] + RSYNC_FLAGS + [
            "--files-from=-", "--from0", "--relative",
            str(src_dir) + "/", str(dst_dir) + "/",
        ]

    return await run_with_retry(cmd_builder, log_path, stdin_bytes=files_bytes)


async def main():
    parser = argparse.ArgumentParser(description="并行同步 Parquet 与 cutflow .out（Python 版，抗抖增强）")
    parser.add_argument("--src-root", default=SRC_ROOT, help="源 Parquet 根目录")
    parser.add_argument("--dst-parq", default=DST_PARQ, help="目标 Parquet 根目录")
    parser.add_argument("--dst-out",  default=DST_OUT,  help="目标 cutflow .out 根目录")
    parser.add_argument("--jobs", type=int, default=3, help="并发任务数（建议 2–3 稳定）")
    parser.add_argument("--logdir", default="./logs_sync", help="日志目录")
    args = parser.parse_args()

    logdir = Path(args.logdir)
    logdir.mkdir(parents=True, exist_ok=True)

    # 预创建目标目录
    for d in SUBDIRS:
        Path("{}/{}".format(args.dst_parq.rstrip("/"), d)).mkdir(parents=True, exist_ok=True)
        Path("{}/{}".format(args.dst_out.rstrip("/"), d)).mkdir(parents=True, exist_ok=True)

    # 任务：先 3 个 Parquet，再 3 个 .out
    task_specs: List[Tuple[str, str]] = [("parquet", d) for d in SUBDIRS] + [("out", d) for d in SUBDIRS]

    sem = asyncio.Semaphore(args.jobs)
    results: List[Tuple[str, str, int]] = []

    async def runner(kind: str, subdir: str):
        async with sem:
            if kind == "parquet":
                rc = await run_parquet_sync(args.src_root, args.dst_parq, subdir, logdir)
            else:
                rc = await run_out_sync(args.src_root, args.dst_out, subdir, logdir)
            results.append((kind, subdir, rc))

    await asyncio.gather(*(runner(k, d) for k, d in task_specs))

    failed = [(k, d, rc) for (k, d, rc) in results if rc != 0]
    if failed:
        print("\n❌ 有任务失败：")
        for k, d, rc in failed:
            print("  - {}:{} -> rc={}".format(k, d, rc))
        print("请查看日志目录：{}".format(logdir))
        sys.exit(1)
    else:
        print("\n✅ 全部同步任务完成。日志在：{}".format(logdir))


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\n^C 中断，已退出。")

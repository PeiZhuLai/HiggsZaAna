ifeq ($(strip $(RecoEgamma/EgammaMLPhotonProducers)),)
ALL_COMMONRULES += src_RecoEgamma_EgammaMLPhotonProducers_src
src_RecoEgamma_EgammaMLPhotonProducers_src_parent := RecoEgamma/EgammaMLPhotonProducers
src_RecoEgamma_EgammaMLPhotonProducers_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_RecoEgamma_EgammaMLPhotonProducers_src,src/RecoEgamma/EgammaMLPhotonProducers/src,LIBRARY))
RecoEgammaEgammaMLPhotonProducers := self/RecoEgamma/EgammaMLPhotonProducers
RecoEgamma/EgammaMLPhotonProducers := RecoEgammaEgammaMLPhotonProducers
RecoEgammaEgammaMLPhotonProducers_files := $(patsubst src/RecoEgamma/EgammaMLPhotonProducers/src/%,%,$(wildcard $(foreach dir,src/RecoEgamma/EgammaMLPhotonProducers/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
RecoEgammaEgammaMLPhotonProducers_BuildFile    := $(WORKINGDIR)/cache/bf/src/RecoEgamma/EgammaMLPhotonProducers/BuildFile
RecoEgammaEgammaMLPhotonProducers_LOC_USE := self   FWCore/Framework FWCore/PluginManager FWCore/ParameterSet PhysicsTools/ONNXRuntime DataFormats/ParticleFlowCandidate DataFormats/ParticleFlowReco Geometry/Records RecoEcal/EgammaCoreTools DataFormats/EcalRecHit Geometry/CaloGeometry Geometry/CaloTopology CondFormats/EcalObjects CondFormats/DataRecord PhysicsTools/UtilAlgos FWCore/ServiceRegistry MLDataFormats/EgammaCandidates PhysicsTools/NanoAOD
RecoEgammaEgammaMLPhotonProducers_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,RecoEgammaEgammaMLPhotonProducers,RecoEgammaEgammaMLPhotonProducers,$(SCRAMSTORENAME_LIB),src/RecoEgamma/EgammaMLPhotonProducers/src))
RecoEgammaEgammaMLPhotonProducers_PACKAGE := self/src/RecoEgamma/EgammaMLPhotonProducers/src
ALL_PRODS += RecoEgammaEgammaMLPhotonProducers
RecoEgammaEgammaMLPhotonProducers_CLASS := LIBRARY
RecoEgamma/EgammaMLPhotonProducers_forbigobj+=RecoEgammaEgammaMLPhotonProducers
RecoEgammaEgammaMLPhotonProducers_INIT_FUNC        += $$(eval $$(call Library,RecoEgammaEgammaMLPhotonProducers,src/RecoEgamma/EgammaMLPhotonProducers/src,src_RecoEgamma_EgammaMLPhotonProducers_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS),edm))
endif

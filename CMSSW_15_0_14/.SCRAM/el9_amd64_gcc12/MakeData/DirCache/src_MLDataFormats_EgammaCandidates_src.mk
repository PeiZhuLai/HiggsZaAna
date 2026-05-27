ifeq ($(strip $(MLDataFormats/EgammaCandidates)),)
ALL_COMMONRULES += src_MLDataFormats_EgammaCandidates_src
src_MLDataFormats_EgammaCandidates_src_parent := MLDataFormats/EgammaCandidates
src_MLDataFormats_EgammaCandidates_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_MLDataFormats_EgammaCandidates_src,src/MLDataFormats/EgammaCandidates/src,LIBRARY))
MLDataFormatsEgammaCandidates := self/MLDataFormats/EgammaCandidates
MLDataFormats/EgammaCandidates := MLDataFormatsEgammaCandidates
MLDataFormatsEgammaCandidates_files := $(patsubst src/MLDataFormats/EgammaCandidates/src/%,%,$(wildcard $(foreach dir,src/MLDataFormats/EgammaCandidates/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
MLDataFormatsEgammaCandidates_BuildFile    := $(WORKINGDIR)/cache/bf/src/MLDataFormats/EgammaCandidates/BuildFile
MLDataFormatsEgammaCandidates_LOC_USE := self   DataFormats/Candidate DataFormats/RecoCandidate FWCore/MessageLogger
MLDataFormatsEgammaCandidates_LCGDICTS  := x 
MLDataFormatsEgammaCandidates_PRE_INIT_FUNC += $$(eval $$(call LCGDict,MLDataFormatsEgammaCandidates,src/MLDataFormats/EgammaCandidates/src/classes.h,src/MLDataFormats/EgammaCandidates/src/classes_def.xml,$(SCRAMSTORENAME_LIB),$(GENREFLEX_ARGS) $(root_EX_FLAGS_GENREFLEX_FAILES_ON_WARNS)))
MLDataFormatsEgammaCandidates_EX_LIB   := MLDataFormatsEgammaCandidates
MLDataFormatsEgammaCandidates_EX_USE   := $(foreach d,$(MLDataFormatsEgammaCandidates_LOC_USE),$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
MLDataFormatsEgammaCandidates_PACKAGE := self/src/MLDataFormats/EgammaCandidates/src
ALL_PRODS += MLDataFormatsEgammaCandidates
MLDataFormatsEgammaCandidates_CLASS := LIBRARY
MLDataFormats/EgammaCandidates_forbigobj+=MLDataFormatsEgammaCandidates
MLDataFormatsEgammaCandidates_INIT_FUNC        += $$(eval $$(call Library,MLDataFormatsEgammaCandidates,src/MLDataFormats/EgammaCandidates/src,src_MLDataFormats_EgammaCandidates_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS),))
endif

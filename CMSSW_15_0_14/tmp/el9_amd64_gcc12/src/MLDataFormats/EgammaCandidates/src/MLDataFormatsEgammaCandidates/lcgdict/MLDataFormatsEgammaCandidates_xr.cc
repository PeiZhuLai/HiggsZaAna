// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdIel9_amd64_gcc12dIsrcdIMLDataFormatsdIEgammaCandidatesdIsrcdIMLDataFormatsEgammaCandidatesdIlcgdictdIMLDataFormatsEgammaCandidates_xr
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "src/MLDataFormats/EgammaCandidates/src/classes.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary();
   static void edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(Long_t size, void *p);
   static void delete_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);
   static void deleteArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);
   static void destruct_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*)
   {
      ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >", ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >::Class_Version(), "DataFormats/Common/interface/Ref.h", 286,
                  typeid(::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >) );
      instance.SetNew(&new_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetNewArray(&newArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDelete(&delete_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDestructor(&destruct_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >","edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >","edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr))->GetClass();
      edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefToBaselErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLRefToBaselErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLRefToBaselErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefToBase<reco::MLPhoton>*)
   {
      ::edm::RefToBase<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefToBase<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefToBase<reco::MLPhoton>", ::edm::RefToBase<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/RefToBase.h", 70,
                  typeid(::edm::RefToBase<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefToBaselErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::RefToBase<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLRefToBaselErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLRefToBaselErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLRefToBaselErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefToBaselErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLRefToBaselErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefToBase<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::RefToBase<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::RefToBase<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefToBaselErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::RefToBase<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLRefToBaselErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefToBaselErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_Dictionary();
   static void edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(Long_t size, void *p);
   static void delete_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p);
   static void deleteArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p);
   static void destruct_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefProd<vector<reco::MLPhoton> >*)
   {
      ::edm::RefProd<vector<reco::MLPhoton> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefProd<vector<reco::MLPhoton> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefProd<vector<reco::MLPhoton> >", ::edm::RefProd<vector<reco::MLPhoton> >::Class_Version(), "DataFormats/Common/interface/RefProd.h", 55,
                  typeid(::edm::RefProd<vector<reco::MLPhoton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::RefProd<vector<reco::MLPhoton> >) );
      instance.SetNew(&new_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetNewArray(&newArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDelete(&delete_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDestructor(&destruct_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::RefProd<vector<reco::MLPhoton> >","edm::RefProd<std::vector<reco::MLPhoton> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::RefProd<vector<reco::MLPhoton> >","edm::RefProd<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefProd<vector<reco::MLPhoton> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::RefProd<vector<reco::MLPhoton> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::RefProd<vector<reco::MLPhoton> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::RefProd<vector<reco::MLPhoton> >*>(nullptr))->GetClass();
      edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary();
   static void edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(Long_t size, void *p);
   static void delete_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);
   static void deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);
   static void destruct_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*)
   {
      ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >", ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >::Class_Version(), "DataFormats/Common/interface/RefVector.h", 32,
                  typeid(::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >) );
      instance.SetNew(&new_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetNewArray(&newArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDelete(&delete_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);
      instance.SetDestructor(&destruct_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >","edm::RefVector<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >","edm::RefVector<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(nullptr))->GetClass();
      edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefToBaseVector<reco::MLPhoton>*)
   {
      ::edm::RefToBaseVector<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefToBaseVector<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefToBaseVector<reco::MLPhoton>", ::edm::RefToBaseVector<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/RefToBaseVector.h", 31,
                  typeid(::edm::RefToBaseVector<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::RefToBaseVector<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefToBaseVector<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::RefToBaseVector<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::RefToBaseVector<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::RefToBaseVector<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefToBaseVectorlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_Dictionary();
   static void edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*)
   {
      ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >", ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >) );
      instance.SetNew(&new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >","edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >","edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> >, unsigned int> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*>(nullptr))->GetClass();
      edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_Dictionary();
   static void edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*)
   {
      ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >", ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >) );
      instance.SetNew(&new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >","edm::Wrapper<edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton>,float,unsigned int> > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >","edm::Wrapper<edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, float, unsigned int> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*>(nullptr))->GetClass();
      edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_Dictionary();
   static void edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*)
   {
      ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >", ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >) );
      instance.SetNew(&new_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >","edm::Wrapper<edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >","edm::Wrapper<edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*>(nullptr))->GetClass();
      edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_Dictionary();
   static void edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*)
   {
      ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >", ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >) );
      instance.SetNew(&new_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*>(nullptr))->GetClass();
      edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::Wrapper<edm::RefVector<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::Wrapper<edm::RefVector<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_Dictionary();
   static void edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(Long_t size, void *p);
   static void delete_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p);
   static void deleteArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p);
   static void destruct_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Wrapper<vector<reco::MLPhoton> >*)
   {
      ::edm::Wrapper<vector<reco::MLPhoton> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Wrapper<vector<reco::MLPhoton> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Wrapper<vector<reco::MLPhoton> >", ::edm::Wrapper<vector<reco::MLPhoton> >::Class_Version(), "DataFormats/Common/interface/Wrapper.h", 25,
                  typeid(::edm::Wrapper<vector<reco::MLPhoton> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Wrapper<vector<reco::MLPhoton> >) );
      instance.SetNew(&new_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetNewArray(&newArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDelete(&delete_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR);
      instance.SetDestructor(&destruct_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<vector<reco::MLPhoton> >","edm::Wrapper<std::vector<reco::MLPhoton> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::Wrapper<vector<reco::MLPhoton> >","edm::Wrapper<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Wrapper<vector<reco::MLPhoton> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Wrapper<vector<reco::MLPhoton> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<vector<reco::MLPhoton> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Wrapper<vector<reco::MLPhoton> >*>(nullptr))->GetClass();
      edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void delete_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::BaseHolder<reco::MLPhoton>*)
   {
      ::edm::reftobase::BaseHolder<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::BaseHolder<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::BaseHolder<reco::MLPhoton>", ::edm::reftobase::BaseHolder<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/BaseHolder.h", 28,
                  typeid(::edm::reftobase::BaseHolder<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::BaseHolder<reco::MLPhoton>) );
      instance.SetDelete(&delete_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::BaseHolder<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::BaseHolder<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::BaseHolder<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::BaseHolder<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/Holder.h", 16,
                  typeid(::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::Holder<reco::Candidate,reco::MLPhotonRef>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::Holder<reco::Candidate, edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::IndirectHolder<reco::MLPhoton>*)
   {
      ::edm::reftobase::IndirectHolder<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::IndirectHolder<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::IndirectHolder<reco::MLPhoton>", ::edm::reftobase::IndirectHolder<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/IndirectHolder.h", 26,
                  typeid(::edm::reftobase::IndirectHolder<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::IndirectHolder<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::IndirectHolder<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::IndirectHolder<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::IndirectHolder<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::IndirectHolder<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/RefHolder_.h", 19,
                  typeid(::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::RefHolder<reco::MLPhotonRef>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::RefHolder<edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void delete_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*)
   {
      ::edm::reftobase::BaseVectorHolder<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::BaseVectorHolder<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::BaseVectorHolder<reco::MLPhoton>", ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/BaseVectorHolder.h", 15,
                  typeid(::edm::reftobase::BaseVectorHolder<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::BaseVectorHolder<reco::MLPhoton>) );
      instance.SetDelete(&delete_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/VectorHolder.h", 14,
                  typeid(::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::VectorHolder<reco::Candidate,reco::MLPhotonRefVector>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::VectorHolder<reco::Candidate, edm::RefVector<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLPtrlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLPtrlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLPtrlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLPtrlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLPtrlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLPtrlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLPtrlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::Ptr<reco::MLPhoton>*)
   {
      ::edm::Ptr<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::Ptr<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::Ptr<reco::MLPhoton>", ::edm::Ptr<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/Ptr.h", 40,
                  typeid(::edm::Ptr<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLPtrlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::Ptr<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLPtrlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLPtrlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLPtrlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLPtrlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLPtrlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::Ptr<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::Ptr<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::Ptr<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLPtrlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::Ptr<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLPtrlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLPtrlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLPtrVectorlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLPtrVectorlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLPtrVectorlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::PtrVector<reco::MLPhoton>*)
   {
      ::edm::PtrVector<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::PtrVector<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::PtrVector<reco::MLPhoton>", ::edm::PtrVector<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/PtrVector.h", 125,
                  typeid(::edm::PtrVector<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLPtrVectorlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::PtrVector<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLPtrVectorlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLPtrVectorlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLPtrVectorlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLPtrVectorlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLPtrVectorlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::PtrVector<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::PtrVector<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::PtrVector<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLPtrVectorlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::PtrVector<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLPtrVectorlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLPtrVectorlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/RefVectorHolder.h", 16,
                  typeid(::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::RefVectorHolder<reco::MLPhotonRefVector>"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::reftobase::RefVectorHolder<edm::RefVector<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLRefToBaseProdlErecocLcLMLPhotongR_Dictionary();
   static void edmcLcLRefToBaseProdlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p);
   static void deleteArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p);
   static void destruct_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::RefToBaseProd<reco::MLPhoton>*)
   {
      ::edm::RefToBaseProd<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::RefToBaseProd<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("edm::RefToBaseProd<reco::MLPhoton>", ::edm::RefToBaseProd<reco::MLPhoton>::Class_Version(), "DataFormats/Common/interface/RefToBaseProd.h", 29,
                  typeid(::edm::RefToBaseProd<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLRefToBaseProdlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::RefToBaseProd<reco::MLPhoton>) );
      instance.SetNew(&new_edmcLcLRefToBaseProdlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_edmcLcLRefToBaseProdlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_edmcLcLRefToBaseProdlErecocLcLMLPhotongR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::RefToBaseProd<reco::MLPhoton>*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::RefToBaseProd<reco::MLPhoton>*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::RefToBaseProd<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLRefToBaseProdlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::RefToBaseProd<reco::MLPhoton>*>(nullptr))->GetClass();
      edmcLcLRefToBaseProdlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLRefToBaseProdlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >::Class_Version(), "DataFormats/Common/interface/ValueMap.h", 107,
                  typeid(::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *recocLcLMLPhoton_Dictionary();
   static void recocLcLMLPhoton_TClassManip(TClass*);
   static void *new_recocLcLMLPhoton(void *p = nullptr);
   static void *newArray_recocLcLMLPhoton(Long_t size, void *p);
   static void delete_recocLcLMLPhoton(void *p);
   static void deleteArray_recocLcLMLPhoton(void *p);
   static void destruct_recocLcLMLPhoton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::reco::MLPhoton*)
   {
      ::reco::MLPhoton *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::reco::MLPhoton));
      static ::ROOT::TGenericClassInfo 
         instance("reco::MLPhoton", 0, "MLDataFormats/EgammaCandidates/interface/MLPhoton.h", 11,
                  typeid(::reco::MLPhoton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &recocLcLMLPhoton_Dictionary, isa_proxy, 12,
                  sizeof(::reco::MLPhoton) );
      instance.SetNew(&new_recocLcLMLPhoton);
      instance.SetNewArray(&newArray_recocLcLMLPhoton);
      instance.SetDelete(&delete_recocLcLMLPhoton);
      instance.SetDeleteArray(&deleteArray_recocLcLMLPhoton);
      instance.SetDestructor(&destruct_recocLcLMLPhoton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::reco::MLPhoton*)
   {
      return GenerateInitInstanceLocal(static_cast<::reco::MLPhoton*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::reco::MLPhoton*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *recocLcLMLPhoton_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::reco::MLPhoton*>(nullptr))->GetClass();
      recocLcLMLPhoton_TClassManip(theClass);
   return theClass;
   }

   static void recocLcLMLPhoton_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_Dictionary();
   static void edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p);
   static void destruct_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*)
   {
      ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >", "DataFormats/Common/interface/AssociationMapHelpers.h", 22,
                  typeid(::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >) );
      instance.SetNew(&new_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >","edm::helpers::KeyVal<edm::RefProd<std::vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >","edm::helpers::KeyVal<edm::RefProd<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> > >, edm::RefProd<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*>(nullptr))->GetClass();
      edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*)
   {
      ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >", "DataFormats/Common/interface/AssociationMapHelpers.h", 37,
                  typeid(::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >) );
      instance.SetNew(&new_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >","edm::helpers::Key<edm::RefProd<std::vector<reco::MLPhoton> > >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >","edm::helpers::Key<edm::RefProd<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> > > >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*>(nullptr))->GetClass();
      edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_Dictionary();
   static void edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p);
   static void deleteArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p);
   static void destruct_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*)
   {
      ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >", ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >::Class_Version(), "DataFormats/Common/interface/AssociationMap.h", 48,
                  typeid(::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >) );
      instance.SetNew(&new_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR);
      instance.SetDelete(&delete_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >","edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >","edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> >, unsigned int> >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*>(nullptr))->GetClass();
      edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_Dictionary();
   static void edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p);
   static void deleteArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p);
   static void destruct_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*)
   {
      ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >", ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >::Class_Version(), "DataFormats/Common/interface/AssociationMap.h", 48,
                  typeid(::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >) );
      instance.SetNew(&new_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR);
      instance.SetDelete(&delete_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR);

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >","edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton>,float,unsigned int> >"));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >","edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, float, unsigned int> >"));
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*)
   {
      return GenerateInitInstanceLocal(static_cast<::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*>(nullptr))->GetClass();
      edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > : new ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >;
   }
   static void *newArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >[nElements] : new ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      delete (static_cast<::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(p));
   }
   static void deleteArray_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      delete [] (static_cast<::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(p));
   }
   static void destruct_edmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      typedef ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBase<reco::MLPhoton> : new ::edm::RefToBase<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLRefToBaselErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBase<reco::MLPhoton>[nElements] : new ::edm::RefToBase<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::RefToBase<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::RefToBase<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLRefToBaselErecocLcLMLPhotongR(void *p) {
      typedef ::edm::RefToBase<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefToBase<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefProd<vector<reco::MLPhoton> > : new ::edm::RefProd<vector<reco::MLPhoton> >;
   }
   static void *newArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefProd<vector<reco::MLPhoton> >[nElements] : new ::edm::RefProd<vector<reco::MLPhoton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete (static_cast<::edm::RefProd<vector<reco::MLPhoton> >*>(p));
   }
   static void deleteArray_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete [] (static_cast<::edm::RefProd<vector<reco::MLPhoton> >*>(p));
   }
   static void destruct_edmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      typedef ::edm::RefProd<vector<reco::MLPhoton> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefProd<vector<reco::MLPhoton> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > : new ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >;
   }
   static void *newArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >[nElements] : new ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      delete (static_cast<::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(p));
   }
   static void deleteArray_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      delete [] (static_cast<::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >*>(p));
   }
   static void destruct_edmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgR(void *p) {
      typedef ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBaseVector<reco::MLPhoton> : new ::edm::RefToBaseVector<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBaseVector<reco::MLPhoton>[nElements] : new ::edm::RefToBaseVector<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::RefToBaseVector<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::RefToBaseVector<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLRefToBaseVectorlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::RefToBaseVector<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefToBaseVector<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > > : new ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >;
   }
   static void *newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >[nElements] : new ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >*>(p));
   }
   static void destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgRsPgR(void *p) {
      typedef ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > > : new ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >;
   }
   static void *newArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >[nElements] : new ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >*>(p));
   }
   static void destruct_edmcLcLWrapperlEedmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgRsPgR(void *p) {
      typedef ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > > : new ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >;
   }
   static void *newArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >[nElements] : new ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >*>(p));
   }
   static void destruct_edmcLcLWrapperlEedmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgRsPgR(void *p) {
      typedef ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> > : new ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >;
   }
   static void *newArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >[nElements] : new ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >*>(p));
   }
   static void destruct_edmcLcLWrapperlEedmcLcLRefToBaseVectorlErecocLcLMLPhotongRsPgR(void *p) {
      typedef ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLWrapperlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<vector<reco::MLPhoton> > : new ::edm::Wrapper<vector<reco::MLPhoton> >;
   }
   static void *newArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Wrapper<vector<reco::MLPhoton> >[nElements] : new ::edm::Wrapper<vector<reco::MLPhoton> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete (static_cast<::edm::Wrapper<vector<reco::MLPhoton> >*>(p));
   }
   static void deleteArray_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      delete [] (static_cast<::edm::Wrapper<vector<reco::MLPhoton> >*>(p));
   }
   static void destruct_edmcLcLWrapperlEvectorlErecocLcLMLPhotongRsPgR(void *p) {
      typedef ::edm::Wrapper<vector<reco::MLPhoton> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Wrapper<vector<reco::MLPhoton> >

namespace ROOT {
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::reftobase::BaseHolder<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::reftobase::BaseHolder<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLBaseHolderlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::reftobase::BaseHolder<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::BaseHolder<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLHolderlErecocLcLCandidatecOedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::IndirectHolder<reco::MLPhoton> : new ::edm::reftobase::IndirectHolder<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::IndirectHolder<reco::MLPhoton>[nElements] : new ::edm::reftobase::IndirectHolder<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::reftobase::IndirectHolder<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::reftobase::IndirectHolder<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLIndirectHolderlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::reftobase::IndirectHolder<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::IndirectHolder<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLRefHolderlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::reftobase::BaseVectorHolder<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLBaseVectorHolderlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::reftobase::BaseVectorHolder<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::BaseVectorHolder<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLVectorHolderlErecocLcLCandidatecOedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLPtrlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Ptr<reco::MLPhoton> : new ::edm::Ptr<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLPtrlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::Ptr<reco::MLPhoton>[nElements] : new ::edm::Ptr<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLPtrlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::Ptr<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLPtrlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::Ptr<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLPtrlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::Ptr<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::Ptr<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::PtrVector<reco::MLPhoton> : new ::edm::PtrVector<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLPtrVectorlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::PtrVector<reco::MLPhoton>[nElements] : new ::edm::PtrVector<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::PtrVector<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::PtrVector<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLPtrVectorlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::PtrVector<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::PtrVector<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLreftobasecLcLRefVectorHolderlEedmcLcLRefVectorlEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBaseProd<reco::MLPhoton> : new ::edm::RefToBaseProd<reco::MLPhoton>;
   }
   static void *newArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::RefToBaseProd<reco::MLPhoton>[nElements] : new ::edm::RefToBaseProd<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<::edm::RefToBaseProd<reco::MLPhoton>*>(p));
   }
   static void deleteArray_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<::edm::RefToBaseProd<reco::MLPhoton>*>(p));
   }
   static void destruct_edmcLcLRefToBaseProdlErecocLcLMLPhotongR(void *p) {
      typedef ::edm::RefToBaseProd<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::RefToBaseProd<reco::MLPhoton>

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLValueMaplEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_recocLcLMLPhoton(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::reco::MLPhoton : new ::reco::MLPhoton;
   }
   static void *newArray_recocLcLMLPhoton(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::reco::MLPhoton[nElements] : new ::reco::MLPhoton[nElements];
   }
   // Wrapper around operator delete
   static void delete_recocLcLMLPhoton(void *p) {
      delete (static_cast<::reco::MLPhoton*>(p));
   }
   static void deleteArray_recocLcLMLPhoton(void *p) {
      delete [] (static_cast<::reco::MLPhoton*>(p));
   }
   static void destruct_recocLcLMLPhoton(void *p) {
      typedef ::reco::MLPhoton current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::reco::MLPhoton

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > > : new ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >;
   }
   static void *newArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >[nElements] : new ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p) {
      delete (static_cast<::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*>(p));
   }
   static void deleteArray_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >*>(p));
   }
   static void destruct_edmcLcLhelperscLcLKeyVallEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRcOedmcLcLRefProdlEedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRsPgRsPgR(void *p) {
      typedef ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > > : new ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >;
   }
   static void *newArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >[nElements] : new ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >*>(p));
   }
   static void destruct_edmcLcLhelperscLcLKeylEedmcLcLRefProdlEvectorlErecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > : new ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >;
   }
   static void *newArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >[nElements] : new ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p) {
      delete (static_cast<::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*>(p));
   }
   static void deleteArray_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p) {
      delete [] (static_cast<::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >*>(p));
   }
   static void destruct_edmcLcLAssociationMaplEedmcLcLOneToOnelEvectorlErecocLcLMLPhotongRcOedmcLcLOwnVectorlErecocLcLCandidatecOedmcLcLClonePolicylErecocLcLCandidategRsPgRcOunsignedsPintgRsPgR(void *p) {
      typedef ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >

namespace ROOT {
   // Wrappers around operator new
   static void *new_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > : new ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >;
   }
   static void *newArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >[nElements] : new ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p) {
      delete (static_cast<::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*>(p));
   }
   static void deleteArray_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p) {
      delete [] (static_cast<::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >*>(p));
   }
   static void destruct_edmcLcLAssociationMaplEedmcLcLOneToValuelEvectorlErecocLcLMLPhotongRcOfloatcOunsignedsPintgRsPgR(void *p) {
      typedef ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >

namespace ROOT {
   static TClass *vectorlErecocLcLMLPhotongR_Dictionary();
   static void vectorlErecocLcLMLPhotongR_TClassManip(TClass*);
   static void *new_vectorlErecocLcLMLPhotongR(void *p = nullptr);
   static void *newArray_vectorlErecocLcLMLPhotongR(Long_t size, void *p);
   static void delete_vectorlErecocLcLMLPhotongR(void *p);
   static void deleteArray_vectorlErecocLcLMLPhotongR(void *p);
   static void destruct_vectorlErecocLcLMLPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<reco::MLPhoton>*)
   {
      vector<reco::MLPhoton> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<reco::MLPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<reco::MLPhoton>", -2, "vector", 423,
                  typeid(vector<reco::MLPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlErecocLcLMLPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<reco::MLPhoton>) );
      instance.SetNew(&new_vectorlErecocLcLMLPhotongR);
      instance.SetNewArray(&newArray_vectorlErecocLcLMLPhotongR);
      instance.SetDelete(&delete_vectorlErecocLcLMLPhotongR);
      instance.SetDeleteArray(&deleteArray_vectorlErecocLcLMLPhotongR);
      instance.SetDestructor(&destruct_vectorlErecocLcLMLPhotongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<reco::MLPhoton> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<reco::MLPhoton>","std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<reco::MLPhoton>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlErecocLcLMLPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<reco::MLPhoton>*>(nullptr))->GetClass();
      vectorlErecocLcLMLPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlErecocLcLMLPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlErecocLcLMLPhotongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<reco::MLPhoton> : new vector<reco::MLPhoton>;
   }
   static void *newArray_vectorlErecocLcLMLPhotongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<reco::MLPhoton>[nElements] : new vector<reco::MLPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlErecocLcLMLPhotongR(void *p) {
      delete (static_cast<vector<reco::MLPhoton>*>(p));
   }
   static void deleteArray_vectorlErecocLcLMLPhotongR(void *p) {
      delete [] (static_cast<vector<reco::MLPhoton>*>(p));
   }
   static void destruct_vectorlErecocLcLMLPhotongR(void *p) {
      typedef vector<reco::MLPhoton> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<reco::MLPhoton>

namespace ROOT {
   static TClass *vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary();
   static void vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void deleteArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);
   static void destruct_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*)
   {
      vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", -2, "vector", 423,
                  typeid(vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >) );
      instance.SetNew(&new_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >","std::vector<edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> >, std::allocator<edm::Ref<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton, edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton, std::allocator<reco::MLPhoton> >, reco::MLPhoton> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(nullptr))->GetClass();
      vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > : new vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >;
   }
   static void *newArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements] : new vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete (static_cast<vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void deleteArray_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      delete [] (static_cast<vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >*>(p));
   }
   static void destruct_vectorlEedmcLcLReflEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotoncOedmcLcLrefhelpercLcLFindUsingAdvancelEvectorlErecocLcLMLPhotongRcOrecocLcLMLPhotongRsPgRsPgR(void *p) {
      typedef vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >

namespace {
  void TriggerDictionaryInitialization_MLDataFormatsEgammaCandidates_xr_Impl() {
    static const char* headers[] = {
"0",
nullptr
    };
    static const char* includePaths[] = {
"src",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_15_0_14/src",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/alpaka/1.2.0-94afb7b501b1959f67f6fd74096328d0/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/pcre/8.43-2d141998cfe5424b8f7aff48035cc2da/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/bz2lib/1.0.6-d065ccd79984efc6d4660f410e4c81de/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gsl/2.6-f7574c606b0ce57ff601d3ca9534cd01/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/libuuid/2.34-27ce4c3579b5b1de2808ea9c4cd8ed29/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/xz/5.2.5-87b9f5597eaeb8b5e9cedb5d183d5089/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/zlib/1.2.13-d217cdbdd8d586e845e05946de2796be/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/eigen/3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e-e265b266d2b30c1bebdd883980d0f9d0/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/eigen/3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e-e265b266d2b30c1bebdd883980d0f9d0/include/eigen3",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fmt/10.2.1-e35fd1db5eb3abc8ac0452e8ee427196/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/md5/1.0.0-5b594b264e04ae51e893b1d69a797ec6/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/OpenBLAS/0.3.27-70a9dd2c9f309171934f13e3003b0540/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/tinyxml2/6.2.0-f99ae2781d074227d47e8a3e7c8ec87e/include",
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.32.13-3ea6c37f14e3f39c0c3ce338c09aec10/include/",
"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MLDataFormatsEgammaCandidates_xr dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace reco{class __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  MLPhoton;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace edm{namespace refhelper{template <typename C, typename T> struct __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefTraits.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  FindUsingAdvance;
}}
namespace edm{namespace refhelper{template <typename REFV> struct __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefTraits.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  FindRefVectorUsingAdvance;
}}
namespace edm{namespace refhelper{template <typename C, typename T> struct __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefTraits.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  FindTrait;
}}
namespace edm{namespace refhelper{template <typename C> struct __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefTraits.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  ValueTrait;
}}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/BaseVectorHolder.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  RefToBase;
}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/Ref.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  RefToBaseVector;
}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/Wrapper.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  Wrapper;
}
namespace edm{namespace reftobase{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/BaseHolder.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  BaseHolder;
}}
namespace reco{class __attribute__((annotate("$clingAutoload$DataFormats/Candidate/interface/Candidate.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  Candidate;}
namespace edm{namespace reftobase{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/IndirectHolder.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  IndirectHolder;
}}
namespace edm{namespace reftobase{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefToBaseVector.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  BaseVectorHolder;
}}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/Ptr.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  Ptr;
}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/HolderToVectorTrait_Ptr_specialization.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  PtrVector;
}
namespace edm{template <typename T> class __attribute__((annotate("$clingAutoload$DataFormats/Common/interface/RefToBaseProd.h")))  __attribute__((annotate("$clingAutoload$MLDataFormats/EgammaCandidates/interface/MLPhoton.h")))  RefToBaseProd;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MLDataFormatsEgammaCandidates_xr dictionary payload"

#ifndef ALPAKA_DEFAULT_HOST_MEMORY_ALIGNMENT
  #define ALPAKA_DEFAULT_HOST_MEMORY_ALIGNMENT 128
#endif
#ifndef ALPAKA_DISABLE_VENDOR_RNG
  #define ALPAKA_DISABLE_VENDOR_RNG 1
#endif
#ifndef CMS_DICT_IMPL
  #define CMS_DICT_IMPL 1
#endif
#ifndef _REENTRANT
  #define _REENTRANT 1
#endif
#ifndef GNUSOURCE
  #define GNUSOURCE 1
#endif
#ifndef __STRICT_ANSI__
  #define __STRICT_ANSI__ 1
#endif
#ifndef CMS_MICRO_ARCH
  #define CMS_MICRO_ARCH "x86-64-v3"
#endif
#ifndef GNU_GCC
  #define GNU_GCC 1
#endif
#ifndef _GNU_SOURCE
  #define _GNU_SOURCE 1
#endif
#ifndef TBB_USE_GLIBCXX_VERSION
  #define TBB_USE_GLIBCXX_VERSION 120301
#endif
#ifndef TBB_SUPPRESS_DEPRECATED_MESSAGES
  #define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#endif
#ifndef TBB_PREVIEW_RESUMABLE_TASKS
  #define TBB_PREVIEW_RESUMABLE_TASKS 1
#endif
#ifndef TBB_PREVIEW_TASK_GROUP_EXTENSIONS
  #define TBB_PREVIEW_TASK_GROUP_EXTENSIONS 1
#endif
#ifndef BOOST_SPIRIT_THREADSAFE
  #define BOOST_SPIRIT_THREADSAFE 1
#endif
#ifndef PHOENIX_THREADSAFE
  #define PHOENIX_THREADSAFE 1
#endif
#ifndef BOOST_MATH_DISABLE_STD_FPCLASSIFY
  #define BOOST_MATH_DISABLE_STD_FPCLASSIFY 1
#endif
#ifndef BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
  #define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX 1
#endif
#ifndef CMSSW_GIT_HASH
  #define CMSSW_GIT_HASH "CMSSW_15_0_14"
#endif
#ifndef PROJECT_NAME
  #define PROJECT_NAME "CMSSW"
#endif
#ifndef PROJECT_VERSION
  #define PROJECT_VERSION "CMSSW_15_0_14"
#endif
#ifndef CMSSW_REFLEX_DICT
  #define CMSSW_REFLEX_DICT 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MLDataFormats/EgammaCandidates/interface/MLPhoton.h"
#include "MLDataFormats/EgammaCandidates/interface/MLPhotonFwd.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Rtypes.h"
#include "Math/Cartesian3D.h"
#include "Math/Polar3D.h"
#include "Math/CylindricalEta3D.h"
#include "Math/PxPyPzE4D.h"
#include <boost/cstdint.hpp>
#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >", payloadCode, "@",
"edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> >", payloadCode, "@",
"edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton>,float,unsigned int> >", payloadCode, "@",
"edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> >", payloadCode, "@",
"edm::Ptr<reco::MLPhoton>", payloadCode, "@",
"edm::PtrVector<reco::MLPhoton>", payloadCode, "@",
"edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> >", payloadCode, "@",
"edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >", payloadCode, "@",
"edm::RefProd<std::vector<reco::MLPhoton> >", payloadCode, "@",
"edm::RefProd<vector<reco::MLPhoton> >", payloadCode, "@",
"edm::RefToBase<reco::MLPhoton>", payloadCode, "@",
"edm::RefToBaseProd<reco::MLPhoton>", payloadCode, "@",
"edm::RefToBaseVector<reco::MLPhoton>", payloadCode, "@",
"edm::RefVector<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> >", payloadCode, "@",
"edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> >", payloadCode, "@",
"edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >", payloadCode, "@",
"edm::Wrapper<edm::AssociationMap<edm::OneToOne<vector<reco::MLPhoton>,edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> >,unsigned int> > >", payloadCode, "@",
"edm::Wrapper<edm::AssociationMap<edm::OneToValue<std::vector<reco::MLPhoton>,float,unsigned int> > >", payloadCode, "@",
"edm::Wrapper<edm::AssociationMap<edm::OneToValue<vector<reco::MLPhoton>,float,unsigned int> > >", payloadCode, "@",
"edm::Wrapper<edm::RefToBaseVector<reco::MLPhoton> >", payloadCode, "@",
"edm::Wrapper<edm::RefVector<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::Wrapper<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::Wrapper<edm::ValueMap<edm::Ref<std::vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<std::vector<reco::MLPhoton>,reco::MLPhoton> > > >", payloadCode, "@",
"edm::Wrapper<edm::ValueMap<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > > >", payloadCode, "@",
"edm::Wrapper<std::vector<reco::MLPhoton> >", payloadCode, "@",
"edm::Wrapper<vector<reco::MLPhoton> >", payloadCode, "@",
"edm::helpers::Key<edm::RefProd<std::vector<reco::MLPhoton> > >", payloadCode, "@",
"edm::helpers::Key<edm::RefProd<vector<reco::MLPhoton> > >", payloadCode, "@",
"edm::helpers::KeyVal<edm::RefProd<std::vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >", payloadCode, "@",
"edm::helpers::KeyVal<edm::RefProd<vector<reco::MLPhoton> >,edm::RefProd<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > >", payloadCode, "@",
"edm::reftobase::BaseHolder<reco::MLPhoton>", payloadCode, "@",
"edm::reftobase::BaseVectorHolder<reco::MLPhoton>", payloadCode, "@",
"edm::reftobase::Holder<reco::Candidate,edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::reftobase::Holder<reco::Candidate,reco::MLPhotonRef>", payloadCode, "@",
"edm::reftobase::IndirectHolder<reco::MLPhoton>", payloadCode, "@",
"edm::reftobase::RefHolder<edm::Ref<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::reftobase::RefHolder<reco::MLPhotonRef>", payloadCode, "@",
"edm::reftobase::RefVectorHolder<edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::reftobase::RefVectorHolder<reco::MLPhotonRefVector>", payloadCode, "@",
"edm::reftobase::VectorHolder<reco::Candidate,edm::RefVector<vector<reco::MLPhoton>,reco::MLPhoton,edm::refhelper::FindUsingAdvance<vector<reco::MLPhoton>,reco::MLPhoton> > >", payloadCode, "@",
"edm::reftobase::VectorHolder<reco::Candidate,reco::MLPhotonRefVector>", payloadCode, "@",
"reco::MLPhoton", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MLDataFormatsEgammaCandidates_xr",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MLDataFormatsEgammaCandidates_xr_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MLDataFormatsEgammaCandidates_xr_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MLDataFormatsEgammaCandidates_xr() {
  TriggerDictionaryInitialization_MLDataFormatsEgammaCandidates_xr_Impl();
}

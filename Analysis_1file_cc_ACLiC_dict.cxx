// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIMartinadIDesktopdICMSdIGIT_HUBdIdisplaced_11dIAnalysis_1file_cc_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "./Analysis_1file.cc"

// Header files passed via #pragma extra_include

   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *ROOT_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("ROOT", 0 /*version*/, "Rtypes.h", 144,
                     ::ROOT::DefineBehavior((void*)0,(void*)0),
                     &ROOT_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *ROOT_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }

namespace ROOT {
   static void *new_Analysis_1file(void *p = 0);
   static void *newArray_Analysis_1file(Long_t size, void *p);
   static void delete_Analysis_1file(void *p);
   static void deleteArray_Analysis_1file(void *p);
   static void destruct_Analysis_1file(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Analysis_1file*)
   {
      ::Analysis_1file *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Analysis_1file >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Analysis_1file", ::Analysis_1file::Class_Version(), "Analysis_1file.h", 38,
                  typeid(::Analysis_1file), DefineBehavior(ptr, ptr),
                  &::Analysis_1file::Dictionary, isa_proxy, 4,
                  sizeof(::Analysis_1file) );
      instance.SetNew(&new_Analysis_1file);
      instance.SetNewArray(&newArray_Analysis_1file);
      instance.SetDelete(&delete_Analysis_1file);
      instance.SetDeleteArray(&deleteArray_Analysis_1file);
      instance.SetDestructor(&destruct_Analysis_1file);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Analysis_1file*)
   {
      return GenerateInitInstanceLocal((::Analysis_1file*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Analysis_1file*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Analysis_1file::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Analysis_1file::Class_Name()
{
   return "Analysis_1file";
}

//______________________________________________________________________________
const char *Analysis_1file::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Analysis_1file*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Analysis_1file::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Analysis_1file*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Analysis_1file::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Analysis_1file*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Analysis_1file::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Analysis_1file*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Analysis_1file::Streamer(TBuffer &R__b)
{
   // Stream an object of class Analysis_1file.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Analysis_1file::Class(),this);
   } else {
      R__b.WriteClassBuffer(Analysis_1file::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Analysis_1file(void *p) {
      return  p ? new(p) ::Analysis_1file : new ::Analysis_1file;
   }
   static void *newArray_Analysis_1file(Long_t nElements, void *p) {
      return p ? new(p) ::Analysis_1file[nElements] : new ::Analysis_1file[nElements];
   }
   // Wrapper around operator delete
   static void delete_Analysis_1file(void *p) {
      delete ((::Analysis_1file*)p);
   }
   static void deleteArray_Analysis_1file(void *p) {
      delete [] ((::Analysis_1file*)p);
   }
   static void destruct_Analysis_1file(void *p) {
      typedef ::Analysis_1file current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Analysis_1file

namespace {
  void TriggerDictionaryInitialization_Analysis_1file_cc_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./Analysis_1file.cc",
0
    };
    static const char* includePaths[] = {
"/Users/Martina/root-6.04.18/include",
"/Users/Martina/root-6.04.18/etc",
"/Users/Martina/root-6.04.18/include",
"/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../include/c++/v1",
"/Users/Martina/root-6.04.18/etc/cling",
"/Users/Martina/root-6.04.18",
"/Users/Martina/root-6.04.18/",
"/Users/Martina/Desktop/CMS/GIT_HUB/displaced_11/",
"graf2d/freetype/src/freetype-2.3.12/include",
"/Users/Martina/root-6.04.18/include",
"/Users/Martina/Desktop/CMS/GIT_HUB/displaced_11/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@/Users/Martina/Desktop/CMS/GIT_HUB/displaced_11/./Analysis_1file.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$./Analysis_1file.cc")))  Analysis_1file;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./Analysis_1file.cc"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"Analysis_1file", payloadCode, "@",
"Analysis_1file::fgIsA", payloadCode, "@",
"ROOT::GenerateInitInstance", payloadCode, "@",
"pigreco", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Analysis_1file_cc_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Analysis_1file_cc_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Analysis_1file_cc_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Analysis_1file_cc_ACLiC_dict() {
  TriggerDictionaryInitialization_Analysis_1file_cc_ACLiC_dict_Impl();
}

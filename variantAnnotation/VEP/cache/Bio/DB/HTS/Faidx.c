/*
 * This file was generated automatically by ExtUtils::ParseXS version 3.28 from the
 * contents of Faidx.xs. Do not edit this file, edit Faidx.xs instead.
 *
 *    ANY CHANGES MADE HERE WILL BE LOST!
 *
 */

#line 1 "lib/Bio/DB/HTS/Faidx.xs"
/*
 * Copyright [2015-2017] EMBL-European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <stdio.h>
#include <string.h>

#include "faidx.h"
#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif


// Code is written to use a blessed int pointer to this strut as an object
// You cannot use Data::Dumper to inspect the Faidx object. Sorry
typedef struct
{
  char* path;
  faidx_t* index;
} Faidx ;


SV* new(const char * classname, const char * path)
{
  Faidx   * faidx;
  SV      * obj;
  SV      * obj_ref;
  faidx_t * fai;

  Newx(faidx, 1, Faidx);

  fai = fai_load(path);
  faidx->path = savepv(path);
  faidx->index = fai;
  obj = newSViv((IV)faidx);
  obj_ref = newRV_noinc(obj);
  sv_bless(obj_ref, gv_stashpv(classname, GV_ADD));
  SvREADONLY_on(obj);

  return obj_ref;
}


void get_sequence(SV* obj, SV* location, SV** seq, int* seq_len)
{
  faidx_t *fai;
  char* char_seq;

  *seq = newSVpvn("",0);
  *seq_len = 0;

  fai = ((Faidx*)SvIV(SvRV(obj)))->index;
  //Fetch sequence
  char_seq = fai_fetch(fai, SvPV(location, PL_na), seq_len);

  //Push into a SV
  sv_catpv(*seq, char_seq);
  sv_2mortal(*seq);
  //Free the buffer created by faidx
  free(char_seq);
}

void get_sequence2(SV* obj, SV* seq_id, int start, int end, SV** seq, int* seq_len)
{
  faidx_t *fai;
  char* char_seq;

  *seq = newSVpvn("",0);
  *seq_len = 0;

  fai = ((Faidx*)SvIV(SvRV(obj)))->index;
  //Fetch sequence
  char_seq = faidx_fetch_seq(fai, SvPV(seq_id, PL_na), start, end, seq_len);

  //Push into a SV
  sv_catpv(*seq, char_seq);
  sv_2mortal(*seq);
  //Free the buffer created by faidx
  free(char_seq);
}


int has_sequence(SV* obj, SV* seq_id)
{
  int has_seq=-1 ;
  has_seq = faidx_has_seq(((Faidx*)SvIV(SvRV(obj)))->index, SvPV(seq_id, PL_na));
  return has_seq;
}


int length(SV* obj, char* seq_id)
{
    int length = 0 ;
    faidx_t *fai = ((Faidx*)SvIV(SvRV(obj)))->index ;
    length = faidx_seq_len(fai, seq_id) ;
    return length ;
}


void DESTROY(SV* obj)
{
  Faidx* faidx = (Faidx*)SvIV(SvRV(obj));
  Safefree(faidx->path);
  fai_destroy(faidx->index);
  Safefree(faidx);
}


#line 135 "lib/Bio/DB/HTS/Faidx.c"
#ifndef PERL_UNUSED_VAR
#  define PERL_UNUSED_VAR(var) if (0) var = var
#endif

#ifndef dVAR
#  define dVAR		dNOOP
#endif


/* This stuff is not part of the API! You have been warned. */
#ifndef PERL_VERSION_DECIMAL
#  define PERL_VERSION_DECIMAL(r,v,s) (r*1000000 + v*1000 + s)
#endif
#ifndef PERL_DECIMAL_VERSION
#  define PERL_DECIMAL_VERSION \
	  PERL_VERSION_DECIMAL(PERL_REVISION,PERL_VERSION,PERL_SUBVERSION)
#endif
#ifndef PERL_VERSION_GE
#  define PERL_VERSION_GE(r,v,s) \
	  (PERL_DECIMAL_VERSION >= PERL_VERSION_DECIMAL(r,v,s))
#endif
#ifndef PERL_VERSION_LE
#  define PERL_VERSION_LE(r,v,s) \
	  (PERL_DECIMAL_VERSION <= PERL_VERSION_DECIMAL(r,v,s))
#endif

/* XS_INTERNAL is the explicit static-linkage variant of the default
 * XS macro.
 *
 * XS_EXTERNAL is the same as XS_INTERNAL except it does not include
 * "STATIC", ie. it exports XSUB symbols. You probably don't want that
 * for anything but the BOOT XSUB.
 *
 * See XSUB.h in core!
 */


/* TODO: This might be compatible further back than 5.10.0. */
#if PERL_VERSION_GE(5, 10, 0) && PERL_VERSION_LE(5, 15, 1)
#  undef XS_EXTERNAL
#  undef XS_INTERNAL
#  if defined(__CYGWIN__) && defined(USE_DYNAMIC_LOADING)
#    define XS_EXTERNAL(name) __declspec(dllexport) XSPROTO(name)
#    define XS_INTERNAL(name) STATIC XSPROTO(name)
#  endif
#  if defined(__SYMBIAN32__)
#    define XS_EXTERNAL(name) EXPORT_C XSPROTO(name)
#    define XS_INTERNAL(name) EXPORT_C STATIC XSPROTO(name)
#  endif
#  ifndef XS_EXTERNAL
#    if defined(HASATTRIBUTE_UNUSED) && !defined(__cplusplus)
#      define XS_EXTERNAL(name) void name(pTHX_ CV* cv __attribute__unused__)
#      define XS_INTERNAL(name) STATIC void name(pTHX_ CV* cv __attribute__unused__)
#    else
#      ifdef __cplusplus
#        define XS_EXTERNAL(name) extern "C" XSPROTO(name)
#        define XS_INTERNAL(name) static XSPROTO(name)
#      else
#        define XS_EXTERNAL(name) XSPROTO(name)
#        define XS_INTERNAL(name) STATIC XSPROTO(name)
#      endif
#    endif
#  endif
#endif

/* perl >= 5.10.0 && perl <= 5.15.1 */


/* The XS_EXTERNAL macro is used for functions that must not be static
 * like the boot XSUB of a module. If perl didn't have an XS_EXTERNAL
 * macro defined, the best we can do is assume XS is the same.
 * Dito for XS_INTERNAL.
 */
#ifndef XS_EXTERNAL
#  define XS_EXTERNAL(name) XS(name)
#endif
#ifndef XS_INTERNAL
#  define XS_INTERNAL(name) XS(name)
#endif

/* Now, finally, after all this mess, we want an ExtUtils::ParseXS
 * internal macro that we're free to redefine for varying linkage due
 * to the EXPORT_XSUB_SYMBOLS XS keyword. This is internal, use
 * XS_EXTERNAL(name) or XS_INTERNAL(name) in your code if you need to!
 */

#undef XS_EUPXS
#if defined(PERL_EUPXS_ALWAYS_EXPORT)
#  define XS_EUPXS(name) XS_EXTERNAL(name)
#else
   /* default to internal */
#  define XS_EUPXS(name) XS_INTERNAL(name)
#endif

#ifndef PERL_ARGS_ASSERT_CROAK_XS_USAGE
#define PERL_ARGS_ASSERT_CROAK_XS_USAGE assert(cv); assert(params)

/* prototype to pass -Wmissing-prototypes */
STATIC void
S_croak_xs_usage(const CV *const cv, const char *const params);

STATIC void
S_croak_xs_usage(const CV *const cv, const char *const params)
{
    const GV *const gv = CvGV(cv);

    PERL_ARGS_ASSERT_CROAK_XS_USAGE;

    if (gv) {
        const char *const gvname = GvNAME(gv);
        const HV *const stash = GvSTASH(gv);
        const char *const hvname = stash ? HvNAME(stash) : NULL;

        if (hvname)
	    Perl_croak_nocontext("Usage: %s::%s(%s)", hvname, gvname, params);
        else
	    Perl_croak_nocontext("Usage: %s(%s)", gvname, params);
    } else {
        /* Pants. I don't think that it should be possible to get here. */
	Perl_croak_nocontext("Usage: CODE(0x%"UVxf")(%s)", PTR2UV(cv), params);
    }
}
#undef  PERL_ARGS_ASSERT_CROAK_XS_USAGE

#define croak_xs_usage        S_croak_xs_usage

#endif

/* NOTE: the prototype of newXSproto() is different in versions of perls,
 * so we define a portable version of newXSproto()
 */
#ifdef newXS_flags
#define newXSproto_portable(name, c_impl, file, proto) newXS_flags(name, c_impl, file, proto, 0)
#else
#define newXSproto_portable(name, c_impl, file, proto) (PL_Sv=(SV*)newXS(name, c_impl, file), sv_setpv(PL_Sv, proto), (CV*)PL_Sv)
#endif /* !defined(newXS_flags) */

#if PERL_VERSION_LE(5, 21, 5)
#  define newXS_deffile(a,b) Perl_newXS(aTHX_ a,b,file)
#else
#  define newXS_deffile(a,b) Perl_newXS_deffile(aTHX_ a,b)
#endif

#line 279 "lib/Bio/DB/HTS/Faidx.c"

XS_EUPXS(XS_Bio__DB__HTS__Faidx_new); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_new)
{
    dVAR; dXSARGS;
    if (items != 2)
       croak_xs_usage(cv,  "classname, path");
    {
	char*	classname = (char *)SvPV_nolen(ST(0))
;
	char*	path = (char *)SvPV_nolen(ST(1))
;
	SV *	RETVAL;

	RETVAL = new(classname, path);
	RETVAL = sv_2mortal(RETVAL);
	ST(0) = RETVAL;
    }
    XSRETURN(1);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence)
{
    dVAR; dXSARGS;
    if (items != 2)
       croak_xs_usage(cv,  "obj, location");
    {
	SV*	obj = ST(0)
;
	SV*	location = ST(1)
;
	SV*	seq;
	int	length;
#line 142 "lib/Bio/DB/HTS/Faidx.xs"
     get_sequence(obj, location, &seq, &length) ;
#line 317 "lib/Bio/DB/HTS/Faidx.c"
	XSprePUSH;	EXTEND(SP,2);
	PUSHs(sv_newmortal());
	ST(0) = seq;
	PUSHs(sv_newmortal());
	sv_setiv(ST(1), (IV)length);
    }
    XSRETURN(2);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence_no_length); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence_no_length)
{
    dVAR; dXSARGS;
    if (items != 2)
       croak_xs_usage(cv,  "obj, location");
    {
	SV*	obj = ST(0)
;
	SV*	location = ST(1)
;
	SV*	seq;
#line 151 "lib/Bio/DB/HTS/Faidx.xs"
  int seq_len=0 ;
  get_sequence(obj, location, &seq, &seq_len) ;
#line 343 "lib/Bio/DB/HTS/Faidx.c"
	XSprePUSH;	EXTEND(SP,1);
	PUSHs(sv_newmortal());
	ST(0) = seq;
    }
    XSRETURN(1);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence2); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence2)
{
    dVAR; dXSARGS;
    if (items != 4)
       croak_xs_usage(cv,  "obj, seq_id, start, end");
    {
	SV*	obj = ST(0)
;
	SV*	seq_id = ST(1)
;
	int	start = (int)SvIV(ST(2))
;
	int	end = (int)SvIV(ST(3))
;
	SV*	seq;
	int	length;

	get_sequence2(obj, seq_id, start, end, &seq, &length);
	XSprePUSH;	EXTEND(SP,2);
	PUSHs(sv_newmortal());
	ST(0) = seq;
	PUSHs(sv_newmortal());
	sv_setiv(ST(1), (IV)length);
    }
    XSRETURN(2);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence2_no_length); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_sequence2_no_length)
{
    dVAR; dXSARGS;
    if (items != 4)
       croak_xs_usage(cv,  "obj, seq_id, start, end");
    {
	SV*	obj = ST(0)
;
	SV*	seq_id = ST(1)
;
	int	start = (int)SvIV(ST(2))
;
	int	end = (int)SvIV(ST(3))
;
	SV*	seq;
#line 173 "lib/Bio/DB/HTS/Faidx.xs"
  int seq_len=0 ;
  get_sequence2(obj, seq_id, start, end, &seq, &seq_len) ;
#line 400 "lib/Bio/DB/HTS/Faidx.c"
	XSprePUSH;	EXTEND(SP,1);
	PUSHs(sv_newmortal());
	ST(0) = seq;
    }
    XSRETURN(1);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_has_sequence); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_has_sequence)
{
    dVAR; dXSARGS;
    if (items != 2)
       croak_xs_usage(cv,  "obj, seq_id");
    {
	SV*	obj = ST(0)
;
	SV*	seq_id = ST(1)
;
	int	RETVAL;
	dXSTARG;

	RETVAL = has_sequence(obj, seq_id);
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_length); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_length)
{
    dVAR; dXSARGS;
    if (items != 2)
       croak_xs_usage(cv,  "obj, seq_id");
    {
	SV*	obj = ST(0)
;
	char*	seq_id = (char *)SvPV_nolen(ST(1))
;
	int	RETVAL;
	dXSTARG;

	RETVAL = length(obj, seq_id);
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_all_sequence_ids); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_get_all_sequence_ids)
{
    dVAR; dXSARGS;
    if (items != 1)
       croak_xs_usage(cv,  "obj");
    PERL_UNUSED_VAR(ax); /* -Wall */
    SP -= items;
    {
	SV*	obj = ST(0)
;
#line 193 "lib/Bio/DB/HTS/Faidx.xs"
   int num_seqs ;
   int i ;
    const char* faidx_name ;
#line 466 "lib/Bio/DB/HTS/Faidx.c"
#line 197 "lib/Bio/DB/HTS/Faidx.xs"
    num_seqs = 0 ;
    faidx_t *fai = ((Faidx*)SvIV(SvRV(obj)))->index ;
    num_seqs = faidx_nseq(fai) ;
    EXTEND(SP,num_seqs);
    for( i=0 ; i<num_seqs ; i++ )
    {
      faidx_name = faidx_iseq(fai,i) ;
      PUSHs(sv_2mortal(newSVpv(faidx_name,0))) ;
    }
#line 477 "lib/Bio/DB/HTS/Faidx.c"
	PUTBACK;
	return;
    }
}


XS_EUPXS(XS_Bio__DB__HTS__Faidx_DESTROY); /* prototype to pass -Wmissing-prototypes */
XS_EUPXS(XS_Bio__DB__HTS__Faidx_DESTROY)
{
    dVAR; dXSARGS;
    if (items != 1)
       croak_xs_usage(cv,  "obj");
    {
	SV*	obj = ST(0)
;

	DESTROY(obj);
    }
    XSRETURN_EMPTY;
}

#ifdef __cplusplus
extern "C"
#endif
XS_EXTERNAL(boot_Bio__DB__HTS__Faidx); /* prototype to pass -Wmissing-prototypes */
XS_EXTERNAL(boot_Bio__DB__HTS__Faidx)
{
#if PERL_VERSION_LE(5, 21, 5)
    dVAR; dXSARGS;
#else
    dVAR; dXSBOOTARGSXSAPIVERCHK;
#endif
#if (PERL_REVISION == 5 && PERL_VERSION < 9)
    char* file = __FILE__;
#else
    const char* file = __FILE__;
#endif

    PERL_UNUSED_VAR(file);

    PERL_UNUSED_VAR(cv); /* -W */
    PERL_UNUSED_VAR(items); /* -W */
#if PERL_VERSION_LE(5, 21, 5)
    XS_VERSION_BOOTCHECK;
#  ifdef XS_APIVERSION_BOOTCHECK
    XS_APIVERSION_BOOTCHECK;
#  endif
#endif

        (void)newXSproto_portable("Bio::DB::HTS::Faidx::new", XS_Bio__DB__HTS__Faidx_new, file, "$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::get_sequence", XS_Bio__DB__HTS__Faidx_get_sequence, file, "$$$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::get_sequence_no_length", XS_Bio__DB__HTS__Faidx_get_sequence_no_length, file, "$$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::get_sequence2", XS_Bio__DB__HTS__Faidx_get_sequence2, file, "$$$$$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::get_sequence2_no_length", XS_Bio__DB__HTS__Faidx_get_sequence2_no_length, file, "$$$$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::has_sequence", XS_Bio__DB__HTS__Faidx_has_sequence, file, "$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::length", XS_Bio__DB__HTS__Faidx_length, file, "$$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::get_all_sequence_ids", XS_Bio__DB__HTS__Faidx_get_all_sequence_ids, file, "$");
        (void)newXSproto_portable("Bio::DB::HTS::Faidx::DESTROY", XS_Bio__DB__HTS__Faidx_DESTROY, file, "$");
#if PERL_VERSION_LE(5, 21, 5)
#  if PERL_VERSION_GE(5, 9, 0)
    if (PL_unitcheckav)
        call_list(PL_scopestack_ix, PL_unitcheckav);
#  endif
    XSRETURN_YES;
#else
    Perl_xs_boot_epilog(aTHX_ ax);
#endif
}


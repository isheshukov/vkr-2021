#ifndef CL_CONFIG_PUBLIC_H
#define CL_CONFIG_PUBLIC_H

#include "cln/host_cpu.h"
#include "cln/version.h"

/* 
 * FIXME: this should not be exposed to user. Or at least it should be
 * renamed to CL_HAVE_LONGLONG or something like that.
 */
/* compiler supports the `long long' type */
#define HAVE_LONGLONG

/**
 * Alignment of a `void*'. CLN needs it to distinguish between pointers
 * and immediate values.
 */
#define cl_word_alignment 8

/* 
 * Numbers in the heap are stored as "digit" (or "limb" in GMP speak)
 * sequences. A digit is an unsigned int with sizeof(void *)*CHAR_BIT bits.
 * It should be 8 or 16 or 32 or 64 bits. If CLN is sitting on top of GMP
 * it should match mp_limb_t
 */
/* #undef GMP_DEMANDS_UINTD_INT */
/* #undef GMP_DEMANDS_UINTD_LONG */
#define GMP_DEMANDS_UINTD_LONG_LONG

#endif /* CL_CONFIG_PUBLIC_H */


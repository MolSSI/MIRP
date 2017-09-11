/*! \file
 *
 * \brief Pragmas for enabling/disabling compiler warnings
 */

#pragma once


#if defined(__ICC) || defined(__INTEL_COMPILER)

    // pragmas for Intel
    #define PRAGMA_WARNING_POP                                 _Pragma("warning(pop)")
    #define PRAGMA_WARNING_PUSH                                _Pragma("warning(push)")
    #define PRAGMA_WARNING_IGNORE_FP_UNDERFLOW                 _Pragma("warning(disable:239)")
    #define PRAGMA_WARNING_IGNORE_FP_EQUALITY                  _Pragma("warning(disable:1572)")
    #define PRAGMA_WARNING_IGNORE_SIGN_CONVERSION              // todo 

#elif defined(__GNUC__) || defined(__GNUG__)

    // pragmas for GCC
    #define PRAGMA_WARNING_PUSH                                _Pragma("GCC diagnostic push")
    #define PRAGMA_WARNING_POP                                 _Pragma("GCC diagnostic pop")
    #define PRAGMA_WARNING_IGNORE_FP_UNDERFLOW                 // Does not exist for gcc?
    #define PRAGMA_WARNING_IGNORE_FP_EQUALITY                  _Pragma("GCC diagnostic ignored \"-Wfloat-equal\"")
    #define PRAGMA_WARNING_IGNORE_SIGN_CONVERSION              _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"")

#endif

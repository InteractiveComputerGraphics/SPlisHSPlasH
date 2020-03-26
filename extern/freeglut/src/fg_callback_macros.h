/*
 * fg_callback_macros.h
 *
 * The freeglut library callback macro file.
 *
 * Copyright (C) 2016 Vincent Simonetti
 * Creation date: Sat Jan 16 2016
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef FREEGLUT_CALLBACK_MACROS_H
#define FREEGLUT_CALLBACK_MACROS_H

/*
 * ----------------------------------------------------------------------------------------------------------------------
 * There are two sets of macros here. One is for executing window callbacks, the others are for setting window callbacks.
 * ----------------------------------------------------------------------------------------------------------------------
 */

/*
 * Compiler define: FG_COMPILER_SUPPORTS_VA_ARGS: if the compiler supports variadic macros
 */

/* What supports variadic macros based off Wikipedia article on it (GCC-like must support C99 or higher to use variadic macros) */
#if (((defined(__GNUC__) && (__GNUC__ >= 3)) || \
      (defined(__clang__))) && \
        (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L))) || \
    (defined(_MSC_VER) && (_MSC_VER >= 1400)) || \
    (defined(__BORLANDC__) && (__BORLANDC__ >= 0x570)) || \
    (defined(__SUNPRO_C) && (__SUNPRO_C >= 0x530))
#define FG_COMPILER_SUPPORTS_VA_ARGS 1
#else
#define FG_COMPILER_SUPPORTS_VA_ARGS 0
#endif

/*
 * --------------------------
 * Executing window callbacks
 * --------------------------
 *
 * Info:
 *
 * This took a while to figure out, so be sure try to understand what is happening so that you can ensure that whatever you
 * change won't break other areas.
 *
 * If you are just adding a new callback/changing it's argument count, just go to the bottom of the file.
 *
 * This whole file exists purely for the sake of preventing the need to implement additional parsing logic for each callback
 * to pass user arguments. Of course, the necessity to support older compilers means that, as seen in the line above, there
 * is still a requirement to add/modify code to handle callbacks. If freeglut ever requires newer compilers (at minimum, ones
 * that support C99 or higher), code can very slowly be removed from this file. Even better would be if the C standard eventually
 * supports something similar to what GCC has implemented or offers an alternative. Another option is if C++ would be "allowed" by
 * project maintaners, as then templates can be used and function overloading. Ironically, the template would probably look worse
 * then the GCC macro, so maybe it's good to stay as is.
 *
 * Onto the different "versions" of macros:
 *
 * The first is for any compiler that supports C99 by default. It requires each callback to have a specific argument count
 * passthrough macro. The only reason there are specific count macros is so that (see paraghraph below) don't need have their own
 * set of callback macros. Ideally, there would only be ZERO and ONE_OR_MORE. This works by having callback-specific macros call a
 * specific handler macro to return user data (ZERO) or return one or more arguments along with userData (ONE_OR_MORE) where, with
 * variadic macros, it just reuses the arguments.
 *
 * The last macro set is for the poor individual who has to use a compiler that doesn't support C99 by default, or may not support
 * it at all. Stuff like MSVC6... It works by having a specific-count macro that "extracts" each argument to have them reused without
 * the parathesis.
 *
 * There is a 3rd macro set that only worked on GCC/Clang, and thus was removed (last seen in revision e9676fc of the GIT mirror.
 * Not sure at this time what the SVN number is.) as it's a non-standard functionality.
 */

/*
 * EXPAND_WCB() is used as:
 *
 *     EXPAND_WCB( cbname )(( arg_list, userData ))
 *
 * ... where {(arg_list)} is the parameter list and userData is user
 * provided data.
 *
 * This will take the arg_list and extend it by one argument, adding
 * the argument "userData" to the end of the list.
 *
 * In order for this to work, each callback must have a define that
 * properly handles the arguments as needed by the callback.
 * This callback is in the format of EXPAND_WCB_SUB_<cbname>.
 * Helper functions exist for zero to five parameters: EXPAND_WCB_ZERO,
 * EXPAND_WCB_ONE, EXPAND_WCB_TWO, EXPAND_WCB_THREE< EXPAND_WCB_FOUR,
 * and EXPAND_WCB_FIVE. Each handle the callback argument counts.
 *
 * An example for the "Entry" callback, where "Entry" is the cbname:
 * typedef void (* FGCBEntry  )( int );
 * typedef void (* FGCBEntryUC)( int, FGCBUserData );
 * #define EXPAND_WCB_SUB_Entry(args) EXPAND_WCB_ONE args
 */

#if FG_COMPILER_SUPPORTS_VA_ARGS

#define EXPAND_WCB_UNPARAN(...) __VA_ARGS__
#define EXPAND_WCB_ONE_OR_MORE(args, userData) ( EXPAND_WCB_UNPARAN args, userData )

#define EXPAND_WCB_ONE(args, userData) EXPAND_WCB_ONE_OR_MORE( args, userData )
#define EXPAND_WCB_TWO(args, userData) EXPAND_WCB_ONE_OR_MORE( args, userData )
#define EXPAND_WCB_THREE(args, userData) EXPAND_WCB_ONE_OR_MORE( args, userData )
#define EXPAND_WCB_FOUR(args, userData) EXPAND_WCB_ONE_OR_MORE( args, userData )
#define EXPAND_WCB_FIVE(args, userData) EXPAND_WCB_ONE_OR_MORE( args, userData )

#else

#define EXPAND_WCB_EXTRACT_ONE_ARGS(arg1) arg1
#define EXPAND_WCB_EXTRACT_TWO_ARGS(arg1, arg2) arg1, arg2
#define EXPAND_WCB_EXTRACT_THREE_ARGS(arg1, arg2, arg3) arg1, arg2, arg3
#define EXPAND_WCB_EXTRACT_FOUR_ARGS(arg1, arg2, arg3, arg4) arg1, arg2, arg3, arg4
#define EXPAND_WCB_EXTRACT_FIVE_ARGS(arg1, arg2, arg3, arg4, arg5) arg1, arg2, arg3, arg4, arg5

#define EXPAND_WCB_ONE(args, userData) (EXPAND_WCB_EXTRACT_ONE_ARGS args, userData)
#define EXPAND_WCB_TWO(args, userData) (EXPAND_WCB_EXTRACT_TWO_ARGS args, userData)
#define EXPAND_WCB_THREE(args, userData) (EXPAND_WCB_EXTRACT_THREE_ARGS args, userData)
#define EXPAND_WCB_FOUR(args, userData) (EXPAND_WCB_EXTRACT_FOUR_ARGS args, userData)
#define EXPAND_WCB_FIVE(args, userData) (EXPAND_WCB_EXTRACT_FIVE_ARGS args, userData)

#endif

#define EXPAND_WCB_ZERO(args, userData) ( userData )

#define EXPAND_WCB(cbname) EXPAND_WCB_SUB_ ## cbname

/*
 * Freeglut callbacks type definitions macros
 *
 * Every time a callback is updated in fg_internal.h is updated, this needs updated
 * if argument counts change, new callbacks are added, or callbacks are removed.
 */

#define EXPAND_WCB_SUB_Display(args) EXPAND_WCB_ZERO args
#define EXPAND_WCB_SUB_Reshape(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_Position(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_Visibility(args) EXPAND_WCB_ONE args
#define EXPAND_WCB_SUB_Keyboard(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_KeyboardUp(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_Special(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_SpecialUp(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_Mouse(args) EXPAND_WCB_FOUR args
#define EXPAND_WCB_SUB_MouseWheel(args) EXPAND_WCB_FOUR args
#define EXPAND_WCB_SUB_Motion(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_Passive(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_Entry(args) EXPAND_WCB_ONE args
#define EXPAND_WCB_SUB_WindowStatus(args) EXPAND_WCB_ONE args
#define EXPAND_WCB_SUB_Joystick(args) EXPAND_WCB_FOUR args
#define EXPAND_WCB_SUB_OverlayDisplay(args) EXPAND_WCB_ZERO args
#define EXPAND_WCB_SUB_SpaceMotion(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_SpaceRotation(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_SpaceButton(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_Dials(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_ButtonBox(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_TabletMotion(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_TabletButton(args) EXPAND_WCB_FOUR args
#define EXPAND_WCB_SUB_Destroy(args) EXPAND_WCB_ZERO args
#define EXPAND_WCB_SUB_MultiEntry(args) EXPAND_WCB_TWO args
#define EXPAND_WCB_SUB_MultiButton(args) EXPAND_WCB_FIVE args
#define EXPAND_WCB_SUB_MultiMotion(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_MultiPassive(args) EXPAND_WCB_THREE args
#define EXPAND_WCB_SUB_InitContext(args) EXPAND_WCB_ZERO args
#define EXPAND_WCB_SUB_AppStatus(args) EXPAND_WCB_ONE args

/*
 * ------------------------
 * Setting window callbacks
 * ------------------------
 *
 * These originally existed in fg_callbacks.c
 */

/*
 * All of the window-specific callbacks setting methods can be generalized to this:
 */
#define SET_CURRENT_WINDOW_CALLBACK(a)                                    \
do                                                                        \
{                                                                         \
    if( fgStructure.CurrentWindow == NULL )                               \
        return;                                                           \
    SET_WCB( ( *( fgStructure.CurrentWindow ) ), a, callback, userData ); \
} while( 0 )

/*
 * Types need to be defined for callbacks. It's not ideal, but it works for this.
 */
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG0(a,b)                              \
static void fgh##a##FuncCallback( FGCBUserData userData )                 \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)();                                                        \
}
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG1(a,b)                              \
static void fgh##a##FuncCallback( int arg1val, FGCBUserData userData )    \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)( arg1val );                                               \
}
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG2(a,b)                              \
static void fgh##a##FuncCallback( int arg1val, int arg2val, FGCBUserData userData ) \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)( arg1val, arg2val );                                      \
}
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG3_USER(a,b,arg1,arg2,arg3)          \
static void fgh##a##FuncCallback( arg1 arg1val, arg2 arg2val, arg3 arg3val, FGCBUserData userData ) \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)( arg1val, arg2val, arg3val );                             \
}
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG3(a,b) IMPLEMENT_CALLBACK_FUNC_CB_ARG3_USER(a,b,int,int,int)
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG4(a,b)                              \
static void fgh##a##FuncCallback( int arg1val, int arg2val, int arg3val, int arg4val, FGCBUserData userData ) \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)( arg1val, arg2val, arg3val, arg4val );                    \
}
#define IMPLEMENT_CALLBACK_FUNC_CB_ARG5(a,b)                              \
static void fgh##a##FuncCallback( int arg1val, int arg2val, int arg3val, int arg4val, int arg5val, FGCBUserData userData ) \
{                                                                         \
    FGCB##b* callback = (FGCB##b*)&userData;                              \
    (*callback)( arg1val, arg2val, arg3val, arg4val, arg5val );           \
}

/*
 * And almost every time the callback setter function can be implemented with these:
 */
#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT_UCALL(a,b)      \
void FGAPIENTRY glut##a##FuncUcall( FGCB##b##UC callback, FGCBUserData userData ) \
{                                                                         \
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glut"#a"FuncUcall" );             \
    SET_CURRENT_WINDOW_CALLBACK( b );                                     \
}
#define IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,b)                      \
void FGAPIENTRY glut##a##Func( FGCB##b callback )                         \
{                                                                         \
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glut"#a"Func" );                  \
    if( callback )                                                        \
    {                                                                     \
        FGCB##b* reference = &callback;                                   \
        glut##a##FuncUcall( fgh##a##FuncCallback, *((FGCBUserData*)reference) ); \
    }                                                                     \
    else                                                                  \
        glut##a##FuncUcall( NULL, NULL );                                 \
}

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,b)            \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT_UCALL(a,b)      \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,b)

/*
 * Combine _glut and _cb macros:
 */
#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG0(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0_2NAME(a,b)             \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG0(a,b)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,b)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG0(a)                               \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG0(a,a)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,a)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG0_2NAME(a,b)                       \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG0(a,b)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,b)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG1(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG1(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG1(a)                               \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG1(a,a)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,a)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG2(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2_2NAME(a,b)             \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG2(a,b)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,b)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG2(a)                               \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG2(a,a)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,a)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG2_2NAME(a,b)                       \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG2(a,b)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,b)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG3(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3_USER(a,arg1,arg2,arg3) \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG3_USER(a,a,arg1,arg2,arg3)           \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG3(a)                               \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG3(a,a)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,a)

#define IMPLEMENT_GLUT_CALLBACK_FUNC_ARG3_2NAME(a,b)                       \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG3(a,b)                               \
        IMPLEMENT_CALLBACK_FUNC_2NAME_GLUT_BASE(a,b)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG4(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG4(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#define IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG5(a)                     \
        IMPLEMENT_CALLBACK_FUNC_CB_ARG5(a,a)                               \
        IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_2NAME_GLUT(a,a)

#endif /* FREEGLUT_CALLBACK_MACROS_H */

/*** END OF FILE ***/

/*
 * fg_callbacks.c
 *
 * The callbacks setting methods.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Fri Dec 3 1999
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

#include <GL/freeglut.h>
#include "fg_internal.h"

/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

/*
 * Global callbacks.
 */

/* Sets the global idle callback */
void FGAPIENTRY glutIdleFuncUcall( FGCBIdleUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutIdleFuncUcall" );
    fgState.IdleCallback = callback;
    fgState.IdleCallbackData = userData;
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG0(Idle)

/* Creates a timer and sets its callback */
void FGAPIENTRY glutTimerFuncUcall( unsigned int timeOut, FGCBTimerUC callback, int timerID, FGCBUserData userData )
{
    SFG_Timer *timer, *node;

    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutTimerFuncUcall" );

    if( (timer = fgState.FreeTimers.Last) )
    {
        fgListRemove( &fgState.FreeTimers, &timer->Node );
    }
    else
    {
        if( ! (timer = malloc(sizeof(SFG_Timer))) )
            fgError( "Fatal error: "
                     "Memory allocation failure in glutTimerFunc()" );
    }

    timer->Callback     = callback;
    timer->CallbackData = userData;
    timer->ID           = timerID;
    timer->TriggerTime  = fgElapsedTime() + timeOut;

    /* Insert such that timers are sorted by end-time */
    for( node = fgState.Timers.First; node; node = node->Node.Next )
    {
        if( node->TriggerTime > timer->TriggerTime )
            break;
    }

    fgListInsert( &fgState.Timers, &node->Node, &timer->Node );
}

IMPLEMENT_CALLBACK_FUNC_CB_ARG1(Timer, Timer)

void FGAPIENTRY glutTimerFunc( unsigned int timeOut, FGCBTimer callback, int timerID )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutTimerFunc" );
    if( callback )
    {
        FGCBTimer* reference = &callback;
        glutTimerFuncUcall( timeOut, fghTimerFuncCallback, timerID, *((FGCBUserData*)reference) );
    }
    else
        glutTimerFuncUcall( timeOut, NULL, timerID, NULL );
}

/* Deprecated version of glutMenuStatusFunc callback setting method */
void FGAPIENTRY glutMenuStateFunc( FGCBMenuState callback )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutMenuStateFunc" );
    fgState.MenuStateCallback = callback;
}

/* Sets the global menu status callback for the current window */
void FGAPIENTRY glutMenuStatusFuncUcall( FGCBMenuStatusUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutMenuStatusFuncUcall" );
    fgState.MenuStatusCallback = callback;
    fgState.MenuStatusCallbackData = userData;
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG3(MenuStatus)

/*
 * Menu specific callbacks.
 */
/* Callback upon menu destruction */
void FGAPIENTRY glutMenuDestroyFuncUcall( FGCBDestroyUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutMenuDestroyFuncUcall" );
    if( fgStructure.CurrentMenu )
    {
        fgStructure.CurrentMenu->Destroy = callback;
        fgStructure.CurrentMenu->DestroyData = userData;
    }
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG0_2NAME(MenuDestroy, Destroy)

/* Implement all these callback setter functions... */
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(Position)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3_USER(Keyboard,unsigned char,int,int)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3_USER(KeyboardUp,unsigned char,int,int)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3(Special)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3(SpecialUp)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG4(Mouse)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG4(MouseWheel)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(Motion)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2_2NAME(PassiveMotion,Passive)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG1(Entry)
/* glutWMCloseFunc is an alias for glutCloseFunc; both set the window's Destroy callback */
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0_2NAME(Close,Destroy)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0_2NAME(WMClose,Destroy)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0(OverlayDisplay)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG1(WindowStatus)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(ButtonBox)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(Dials)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(TabletMotion)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG4(TabletButton)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG2(MultiEntry)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG5(MultiButton)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3(MultiMotion)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG3(MultiPassive)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG0(InitContext)
IMPLEMENT_CURRENT_WINDOW_CALLBACK_FUNC_ARG1(AppStatus)

/*
 * Sets the Display callback for the current window
 */
void FGAPIENTRY glutDisplayFuncUcall( FGCBDisplayUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutDisplayFuncUcall" );
    if( !callback )
        fgError( "Fatal error in program.  NULL display callback not "
                 "permitted in GLUT 3.0+ or freeglut 2.0.1+" );
    SET_CURRENT_WINDOW_CALLBACK( Display );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG0(Display)

void fghDefaultReshape( int width, int height, FGCBUserData userData )
{
    glViewport( 0, 0, width, height );
}

void FGAPIENTRY glutReshapeFuncUcall( FGCBReshapeUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutReshapeFuncUcall" );
    
    if( !callback )
    {
        callback = fghDefaultReshape;
        userData = NULL;
    }

    SET_CURRENT_WINDOW_CALLBACK( Reshape );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG2(Reshape)

/*
 * Sets the Visibility callback for the current window.
 * NB: the Visibility func is deprecated in favor of the WindowStatus func,
 * which provides more detail. The visibility func callback is implemented
 * as a translation step from the windowStatus func. When the user sets the
 * windowStatus func, any visibility func is overwritten.
 * DEVELOPER NOTE: in the library, only invoke the window status func, this
 * gets automatically translated to the visibility func if thats what the
 * user has set.
 * window status is kind of anemic on win32 as there are no window messages
 * to notify us that the window is covered by other windows or not.
 * Should one want to query this, see
 * http://stackoverflow.com/questions/5445889/get-which-process-window-is-actually-visible-in-c-sharp
 * for an implementation outline (but it would be polling based, not push based).
 */
static void fghVisibility( int status, FGCBUserData userData )
{
    int vis_status;

    FREEGLUT_INTERNAL_ERROR_EXIT_IF_NOT_INITIALISED ( "Visibility Callback" );
    freeglut_return_if_fail( fgStructure.CurrentWindow );

    /* Translate window status func states to visibility states */
    if( ( status == GLUT_HIDDEN)  || ( status == GLUT_FULLY_COVERED) )
        vis_status = GLUT_NOT_VISIBLE;
    else    /* GLUT_FULLY_RETAINED, GLUT_PARTIALLY_RETAINED */
        vis_status = GLUT_VISIBLE;

    INVOKE_WCB( *( fgStructure.CurrentWindow ), Visibility, ( vis_status ) );
}

void FGAPIENTRY glutVisibilityFuncUcall( FGCBVisibilityUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutVisibilityFuncUcall" );

    if ( !callback )
    {
        userData = NULL;
    }

    SET_CURRENT_WINDOW_CALLBACK( Visibility );

    if( callback )
        glutWindowStatusFuncUcall( fghVisibility, NULL );
    else
        glutWindowStatusFuncUcall( NULL, NULL );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG1(Visibility)

/*
 * Sets the joystick callback and polling rate for the current window
 */
void FGAPIENTRY glutJoystickFuncUcall( FGCBJoystickUC callback, int pollInterval, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutJoystickFuncUcall" );
    fgInitialiseJoysticks ();

    if ( (
           fgStructure.CurrentWindow->State.JoystickPollRate <= 0 ||        /* Joystick callback was disabled */
           !FETCH_WCB(*fgStructure.CurrentWindow,Joystick)
         ) &&
         ( 
           callback && ( pollInterval > 0 )                                 /* but is now enabled */
         ) )
        ++fgState.NumActiveJoysticks;
    else if ( ( 
                fgStructure.CurrentWindow->State.JoystickPollRate > 0 &&    /* Joystick callback was enabled */
                FETCH_WCB(*fgStructure.CurrentWindow,Joystick)
              ) &&  
              ( 
                !callback || ( pollInterval <= 0 )                          /* but is now disabled */
              ) )
        --fgState.NumActiveJoysticks;

    SET_CURRENT_WINDOW_CALLBACK( Joystick );
    fgStructure.CurrentWindow->State.JoystickPollRate = pollInterval;

    /* set last poll time such that joystick will be polled asap */
    fgStructure.CurrentWindow->State.JoystickLastPoll = fgElapsedTime();
    if (fgStructure.CurrentWindow->State.JoystickLastPoll < pollInterval)
        fgStructure.CurrentWindow->State.JoystickLastPoll = 0;
    else
        fgStructure.CurrentWindow->State.JoystickLastPoll -= pollInterval;
}

static void fghJoystickFuncCallback( unsigned int buttons, int axis0, int axis1, int axis2, FGCBUserData userData )
{
    FGCBJoystick* callback = (FGCBJoystick*)&userData;
    (*callback)( buttons, axis0, axis1, axis2 );
}

void FGAPIENTRY glutJoystickFunc( FGCBJoystick callback, int pollInterval )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutJoystickFunc" );
    if( callback )
    {
        FGCBJoystick* reference = &callback;
        glutJoystickFuncUcall( fghJoystickFuncCallback, pollInterval, *((FGCBUserData*)reference) );
    }
    else
        glutJoystickFuncUcall( NULL, pollInterval, NULL );
}

/*
 * Sets the spaceball motion callback for the current window
 */
void FGAPIENTRY glutSpaceballMotionFuncUcall( FGCBSpaceMotionUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutSpaceballMotionFuncUcall" );
    fgInitialiseSpaceball();

    SET_CURRENT_WINDOW_CALLBACK( SpaceMotion );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG3_2NAME(SpaceballMotion, SpaceMotion)

/*
 * Sets the spaceball rotate callback for the current window
 */
void FGAPIENTRY glutSpaceballRotateFuncUcall( FGCBSpaceRotationUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutSpaceballRotateFuncUcall" );
    fgInitialiseSpaceball();

    SET_CURRENT_WINDOW_CALLBACK( SpaceRotation );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG3_2NAME(SpaceballRotate, SpaceRotation)

/*
 * Sets the spaceball button callback for the current window
 */
void FGAPIENTRY glutSpaceballButtonFuncUcall( FGCBSpaceButtonUC callback, FGCBUserData userData )
{
    FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutSpaceballButtonFuncUcall" );
    fgInitialiseSpaceball();

    SET_CURRENT_WINDOW_CALLBACK( SpaceButton );
}

IMPLEMENT_GLUT_CALLBACK_FUNC_ARG2_2NAME(SpaceballButton, SpaceButton)

/*** END OF FILE ***/

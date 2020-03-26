/*
 * fg_spaceball_mswin.c
 *
 * Spaceball support for Windows
 *
 * Copyright (c) 2012 Stephen J. Baker. All Rights Reserved.
 * Written by Evan Felix <karcaw at gmail.com>
 * Creation date: Sat Feb 4, 2012
 *
 * Copyright (c) 2014 Jinrong Xie. All Rights Reserved.
 * Written by Jinrong Xie <stonexjr at gmail.com>
 * Modification date: Wed Dec 24, 2014
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

/*
 * Modified by Jinrong Xie <stonexjr at gmail.com> 12/24/2014
 * for Space Navigator support on Windows.
 * This code is enhanced by at least supporting 3Dconnexion's
 * six degree of freedom navigator.
 */

//#if(_WIN32_WINNT >= 0x0501)

#include <GL/freeglut.h>
#include <stdlib.h>
#include "../fg_internal.h"

enum {
    SPNAV_EVENT_ANY,
    SPNAV_EVENT_MOTION_TRANSLATION,
    SPNAV_EVENT_MOTION_ROTATION,
    SPNAV_EVENT_BUTTON  /* includes both press and release */
};

extern int fg_sball_initialized;
unsigned int __fgSpaceKeystate = 0;
RAWINPUTDEVICE __fgSpaceball = { 0x01, 0x08, 0x00, 0x00 };

void fgPlatformInitializeSpaceball(void)
{
    HWND hwnd;
    fg_sball_initialized = 1;
    if (!fgStructure.CurrentWindow)
    {
        fg_sball_initialized = 0;
        return;
    }
    hwnd = fgStructure.CurrentWindow->Window.Handle;

    {
        BOOL ok;
        UINT cbSize = sizeof(__fgSpaceball);
        __fgSpaceball.hwndTarget = hwnd;
        ok = RegisterRawInputDevices(&__fgSpaceball, 1, cbSize);

        if (!ok){
            __fgSpaceball.hwndTarget = NULL;
            fg_sball_initialized = 0;
        }
    }
}

void fgPlatformSpaceballClose(void)
{
    return;
}

int fgPlatformHasSpaceball(void)
{
    return __fgSpaceball.hwndTarget ? 1 : 0;
}

int fgPlatformSpaceballNumButtons(void)
{
    return 0;
}

void fgPlatformSpaceballSetWindow(SFG_Window *window)
{
    return;
}

int fgIsSpaceballWinEvent(HWND hwnd, WPARAM wParam, LPARAM lParam)
{
    return 0;
}

void fgSpaceballHandleWinEvent(HWND hwnd, WPARAM wParam, LPARAM lParam)
{
    #define LOGITECH_VENDOR_ID 0x46d
    HRAWINPUT hRawInput = (HRAWINPUT)lParam;
    UINT inputCode = (UINT)wParam;
    UINT size;
    BYTE *rawInputBuffer;
    PRAWINPUT pRawInput;
    UINT res;
    RID_DEVICE_INFO sRidDeviceInfo;

    if (!fg_sball_initialized)
    {
        fgPlatformInitializeSpaceball();
        if (!fg_sball_initialized)
        {
            return;
        }
    }

    res = GetRawInputData(hRawInput, RID_INPUT, NULL, &size, sizeof(RAWINPUTHEADER));
    if (res == -1)
        return;

    rawInputBuffer = malloc(size * sizeof *rawInputBuffer);
    pRawInput = (PRAWINPUT)rawInputBuffer;

    res = GetRawInputData(hRawInput, RID_INPUT, pRawInput, &size, sizeof(RAWINPUTHEADER));
    if (res == -1)
        return;
    if (pRawInput->header.dwType != RIM_TYPEHID)
        return;

    sRidDeviceInfo.cbSize = sizeof(RID_DEVICE_INFO);
    size = sizeof(RID_DEVICE_INFO);
    res = GetRawInputDeviceInfo(pRawInput->header.hDevice, RIDI_DEVICEINFO, &sRidDeviceInfo, &size);
    if (res == -1)
        return;
    {
        SFG_Window* window = fgWindowByHandle(hwnd);
        if ((window == NULL))
            return;

        if (sRidDeviceInfo.hid.dwVendorId == LOGITECH_VENDOR_ID)
        {
            // Motion data comes in two parts: motion type and
            // displacement/rotation along three axis.
            // Orientation is a right handed coordinate system with
            // X goes right, Y goes up and Z goes towards viewer, e.g.
            // the one used in OpenGL
            if (pRawInput->data.hid.bRawData[0] ==
                SPNAV_EVENT_MOTION_TRANSLATION)
            { // Translation vector
                short* pnData = (short*)(&pRawInput->data.hid.bRawData[1]);
                short X = pnData[0];
                short Y = -pnData[2];
                short Z = pnData[1];
                INVOKE_WCB(*window, SpaceMotion, (X, Y, Z));
            }
            else if (pRawInput->data.hid.bRawData[0] ==
                SPNAV_EVENT_MOTION_ROTATION)
            { // Axis aligned rotation vector
                short* pnData = (short*)(&pRawInput->data.hid.bRawData[1]);
                short rX = pnData[0];
                short rY = -pnData[2];
                short rZ = pnData[1];
                INVOKE_WCB(*window, SpaceRotation, (rX, rY, rZ));
            }
            else if (pRawInput->data.hid.bRawData[0] ==
                SPNAV_EVENT_BUTTON)
            { // State of the keys
                unsigned long dwKeystate = *(unsigned long*)(&pRawInput->data.hid.bRawData[1]);
                unsigned int state = GLUT_UP;
                if (FETCH_WCB(*window, SpaceButton))
                {
                    int i;
                    for (i = 0; i < 32; i++)
                    {
                        unsigned long stateBefore = __fgSpaceKeystate&(1 << i);
                        unsigned long stateNow = dwKeystate&(1 << i);

                        if (stateBefore && !stateNow)
                            INVOKE_WCB(*window, SpaceButton, (stateBefore, GLUT_UP));
                        if (!stateBefore && stateNow)
                            INVOKE_WCB(*window, SpaceButton, (stateNow, GLUT_DOWN));

                    }
                }
                __fgSpaceKeystate = dwKeystate;
            }
        }
    }
}

//#endif

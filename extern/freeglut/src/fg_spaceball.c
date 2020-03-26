/* Spaceball support for Linux.
 * Written by John Tsiombikas <nuclear@member.fsf.org>
 *
 * This code supports 3Dconnexion's 6-dof space-whatever devices.
 * It can communicate with either the proprietary 3Dconnexion daemon (3dxsrv)
 * free spacenavd (http://spacenav.sourceforge.net), through the "standard"
 * magellan X-based protocol.
 */


#include <GL/freeglut.h>
#include "fg_internal.h"

#if( !_WIN32 || _WIN32_WINNT >= 0x0501)

/* -- PRIVATE FUNCTIONS --------------------------------------------------- */

extern void fgPlatformInitializeSpaceball(void);
extern void fgPlatformSpaceballClose(void);
extern int fgPlatformHasSpaceball(void);
extern int fgPlatformSpaceballNumButtons(void);
extern void fgPlatformSpaceballSetWindow(SFG_Window *window);


int fg_sball_initialized = 0;

void fgInitialiseSpaceball(void)
{
    if(fg_sball_initialized != 0) {
        return;
    }

    fgPlatformInitializeSpaceball();
}

void fgSpaceballClose(void)
{
    fgPlatformSpaceballClose();
}

int fgHasSpaceball(void)
{
    if(fg_sball_initialized == 0) {
        fgInitialiseSpaceball();
        if(fg_sball_initialized != 1) {
            fgWarning("fgInitialiseSpaceball failed\n");
            return 0;
        }
    }

    return fgPlatformHasSpaceball();
}

int fgSpaceballNumButtons(void)
{
    if(fg_sball_initialized == 0) {
        fgInitialiseSpaceball();
        if(fg_sball_initialized != 1) {
            fgWarning("fgInitialiseSpaceball failed\n");
            return 0;
        }
    }

    return fgPlatformSpaceballNumButtons();
}

void fgSpaceballSetWindow(SFG_Window *window)
{
    if(fg_sball_initialized == 0) {
        fgInitialiseSpaceball();
        if(fg_sball_initialized != 1) {
            return;
        }
    }

    fgPlatformSpaceballSetWindow(window);
}

#else

void fgInitialiseSpaceball(void)
{
}

void fgSpaceballClose(void)
{
}

int fgHasSpaceball(void)
{
    return 0;
}

int fgSpaceballNumButtons(void)
{
    return 0;
}

void fgSpaceballSetWindow(SFG_Window *window)
{
}

#endif

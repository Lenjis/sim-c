#include <stdio.h>
#include <windows.h>
#include <mmsystem.h>

// Link the winmm.lib library for multimedia functions
#pragma comment(lib, "winmm.lib")

// Timer callback function that gets called when the timer fires
VOID CALLBACK TimerProc(HWND hwnd,  // Handle to the window (not used here)
                        UINT uMsg,  // Message identifier (not used here)
                        UINT_PTR idEvent,  // Timer identifier
                        DWORD dwTime  // Current system time in milliseconds
) {
    // Print the time at which the timer fired
    printf("Timer fired at %d ms.\n", dwTime);
}

int main() {
    UINT uDelay = 1000;    // Timer interval in milliseconds (1 second)
    UINT uResolution = 1;  // Timer resolution in milliseconds
    // Create a periodic timer
    UINT timerID =
        timeSetEvent(uDelay, uResolution, TimerProc, 0, TIME_PERIODIC);
    if (timerID == NULL)  // Check if the timer was created successfully
    {
        printf("Failed to create timer.\n");
        return 1;  // Exit with error code
    }

    Sleep(10000);  // Keep the program running for 10 seconds to allow the timer
                   // to fire
    timeKillEvent(timerID);  // Stop the timer
    return 0;                // Exit successfully
}

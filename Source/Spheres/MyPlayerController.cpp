// MyPlayerController.cpp
#include "MyPlayerController.h"

void AMyPlayerController::BeginPlay()
{
    Super::BeginPlay();

    bEnableTouchEvents = true;
    bEnableClickEvents = true;
    bEnableMouseOverEvents = true;
}

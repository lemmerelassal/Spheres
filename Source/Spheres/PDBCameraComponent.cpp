#include "PDBCameraComponent.h"
#include "Camera/CameraComponent.h"
#include "Components/SceneComponent.h"
#include "GameFramework/PlayerController.h"
#include "Components/InputComponent.h"
#include "Framework/Application/SlateApplication.h"

APDBCameraComponent::APDBCameraComponent()
{
    PrimaryActorTick.bCanEverTick = true;

    // Create pivot point at center
    PivotPoint = CreateDefaultSubobject<USceneComponent>(TEXT("PivotPoint"));
    RootComponent = PivotPoint;

    // Create camera
    Camera = CreateDefaultSubobject<UCameraComponent>(TEXT("Camera"));
    Camera->SetupAttachment(PivotPoint);

    // Default camera settings
    OrbitDistance = 2000.0f;
    RotationSpeed = 1.0f;
    PanSpeed = 5.0f;
    ZoomSpeed = 100.0f;
    MinZoomDistance = 200.0f;
    MaxZoomDistance = 10000.0f;

    CurrentYaw = 0.0f;
    CurrentPitch = -30.0f;
    bLeftMouseWasDown = false;
    bRightMouseWasDown = false;
    bIsFirstLeftFrame = false;
    bIsFirstRightFrame = false;

    // Set initial camera position
    Camera->SetRelativeLocation(FVector(-OrbitDistance, 0, 0));

    // Enable input
    AutoPossessPlayer = EAutoReceiveInput::Player0;
}

void APDBCameraComponent::BeginPlay()
{
    Super::BeginPlay();

    // Get player controller and set view target
    APlayerController* PC = GetWorld()->GetFirstPlayerController();
    if (PC)
    {
        PC->SetViewTarget(this);
        PC->bShowMouseCursor = true;
        
        // Disable touch interface (removes on-screen joysticks)
        PC->ActivateTouchInterface(nullptr);
        
        // Use Game and UI mode
        FInputModeGameAndUI InputMode;
        InputMode.SetHideCursorDuringCapture(false);
        InputMode.SetLockMouseToViewportBehavior(EMouseLockMode::DoNotLock);
        PC->SetInputMode(InputMode);
    }

    // Apply initial rotation
    PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
    
    LastMousePosition = FVector2D::ZeroVector;
    MousePositionOnPress = FVector2D::ZeroVector;
}

void APDBCameraComponent::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);
    
    APlayerController* PC = GetWorld()->GetFirstPlayerController();
    if (!PC) return;

    // Get mouse position from Slate (this is more reliable)
    FVector2D CurrentMousePos = FSlateApplication::Get().GetCursorPos();
    
    // Check mouse buttons - SWAPPED: Right=Rotate, Middle=Pan
    bool bRotateMouseDown = PC->IsInputKeyDown(EKeys::RightMouseButton);
    bool bPanMouseDown = PC->IsInputKeyDown(EKeys::MiddleMouseButton);

    // Handle rotate mouse button state changes (RIGHT MOUSE)
    if (bRotateMouseDown && !bLeftMouseWasDown)
    {
        // Just pressed
        MousePositionOnPress = CurrentMousePos;
        bIsFirstLeftFrame = true;
    }
    else if (!bRotateMouseDown && bLeftMouseWasDown)
    {
        // Just released
        bIsFirstLeftFrame = false;
    }

    // Handle pan mouse button state changes (MIDDLE MOUSE)
    if (bPanMouseDown && !bRightMouseWasDown)
    {
        // Just pressed
        MousePositionOnPress = CurrentMousePos;
        bIsFirstRightFrame = true;
    }
    else if (!bPanMouseDown && bRightMouseWasDown)
    {
        // Just released
        bIsFirstRightFrame = false;
    }

    // Calculate delta
    FVector2D MouseDelta = CurrentMousePos - LastMousePosition;

    // Rotation with RIGHT mouse (skip first frame after press to avoid jump)
    if (bRotateMouseDown && !bIsFirstLeftFrame && MouseDelta.SizeSquared() > 0.01f)
    {
        CurrentYaw += MouseDelta.X * RotationSpeed * 0.5f;
        CurrentPitch = FMath::Clamp(CurrentPitch - MouseDelta.Y * RotationSpeed * 0.5f, -89.0f, 89.0f);
        PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
    }

    // Pan with MIDDLE mouse (skip first frame after press to avoid jump)
    if (bPanMouseDown && !bIsFirstRightFrame && MouseDelta.SizeSquared() > 0.01f)
    {
        FVector RightVector = Camera->GetRightVector();
        FVector UpVector = Camera->GetUpVector();
        FVector NewLocation = GetActorLocation() - RightVector * MouseDelta.X * PanSpeed * 0.5f + UpVector * MouseDelta.Y * PanSpeed * 0.5f;
        SetActorLocation(NewLocation);
    }

    // Clear first frame flags after first frame
    if (bIsFirstLeftFrame)
        bIsFirstLeftFrame = false;
    if (bIsFirstRightFrame)
        bIsFirstRightFrame = false;

    // Update tracking variables
    LastMousePosition = CurrentMousePos;
    bLeftMouseWasDown = bRotateMouseDown;
    bRightMouseWasDown = bPanMouseDown;

    // Mouse wheel zoom
    float WheelAxis = PC->GetInputAnalogKeyState(EKeys::MouseWheelAxis);
    if (FMath::Abs(WheelAxis) > 0.01f)
    {
        OrbitDistance = FMath::Clamp(OrbitDistance - WheelAxis * ZoomSpeed, MinZoomDistance, MaxZoomDistance);
        Camera->SetRelativeLocation(FVector(-OrbitDistance, 0, 0));
    }
}

void APDBCameraComponent::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
    Super::SetupPlayerInputComponent(PlayerInputComponent);
}

void APDBCameraComponent::SetTargetActor(AActor* Target)
{
    if (Target)
    {
        SetActorLocation(Target->GetActorLocation());
    }
}
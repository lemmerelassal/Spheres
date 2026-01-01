#include "PDBCameraComponent.h"
#include "Camera/CameraComponent.h"
#include "Components/SceneComponent.h"
#include "GameFramework/PlayerController.h"
#include "Components/InputComponent.h"
#include "Framework/Application/SlateApplication.h"
#include "Kismet/GameplayStatics.h"
#include "Components/PrimitiveComponent.h"

APDBCameraComponent::APDBCameraComponent()
{
    PrimaryActorTick.bCanEverTick = true;

    // Create pivot point at center
    PivotPoint = CreateDefaultSubobject<USceneComponent>(TEXT("PivotPoint"));
    RootComponent = PivotPoint;

    // Create camera
    Camera = CreateDefaultSubobject<UCameraComponent>(TEXT("Camera"));
    Camera->SetupAttachment(PivotPoint);
    
    // Set camera clipping planes for very large scenes
    Camera->SetConstraintAspectRatio(false);
    Camera->PostProcessSettings.bOverride_DepthOfFieldFstop = false;

    // Default camera settings
    OrbitDistance = 2000.0f;
    RotationSpeed = 1.0f;
    PanSpeed = 5.0f;
    ZoomSpeed = 100.0f;
    MinZoomDistance = 200.0f;
    MaxZoomDistance = 50000000.0f; // Very large scenes supported

    CurrentYaw = 0.0f;
    CurrentPitch = -30.0f;
    bLeftMouseWasDown = false;
    bRightMouseWasDown = false;
    bIsFirstLeftFrame = false;
    bIsFirstRightFrame = false;

    // Initial camera offset
    Camera->SetRelativeLocation(FVector(-OrbitDistance, 0, 0));

    // Enable player input
    AutoPossessPlayer = EAutoReceiveInput::Player0;
}

void APDBCameraComponent::BeginPlay()
{
    Super::BeginPlay();

    APlayerController* PC = GetWorld()->GetFirstPlayerController();
    if (PC)
    {
        PC->SetViewTarget(this);
        PC->bShowMouseCursor = true;
        PC->ActivateTouchInterface(nullptr);

        FInputModeGameAndUI InputMode;
        InputMode.SetHideCursorDuringCapture(false);
        InputMode.SetLockMouseToViewportBehavior(EMouseLockMode::DoNotLock);
        PC->SetInputMode(InputMode);
    }

    PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));

    LastMousePosition = FVector2D::ZeroVector;
    MousePositionOnPress = FVector2D::ZeroVector;
}

void APDBCameraComponent::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    APlayerController* PC = GetWorld()->GetFirstPlayerController();
    if (!PC) return;

    FVector2D CurrentMousePos = FSlateApplication::Get().GetCursorPos();

    // Input: right = rotate, middle = pan
    bool bRotateMouseDown = PC->IsInputKeyDown(EKeys::RightMouseButton);
    bool bPanMouseDown = PC->IsInputKeyDown(EKeys::MiddleMouseButton);

    // Handle right mouse press/release
    if (bRotateMouseDown && !bLeftMouseWasDown)
    {
        MousePositionOnPress = CurrentMousePos;
        bIsFirstLeftFrame = true;
    }
    else if (!bRotateMouseDown && bLeftMouseWasDown)
    {
        bIsFirstLeftFrame = false;
    }

    // Handle middle mouse press/release
    if (bPanMouseDown && !bRightMouseWasDown)
    {
        MousePositionOnPress = CurrentMousePos;
        bIsFirstRightFrame = true;
    }
    else if (!bPanMouseDown && bRightMouseWasDown)
    {
        bIsFirstRightFrame = false;
    }

    FVector2D MouseDelta = CurrentMousePos - LastMousePosition;

    // Rotate
    if (bRotateMouseDown && !bIsFirstLeftFrame && MouseDelta.SizeSquared() > 0.01f)
    {
        CurrentYaw += MouseDelta.X * RotationSpeed * 0.5f;
        CurrentPitch = FMath::Clamp(CurrentPitch - MouseDelta.Y * RotationSpeed * 0.5f, -89.0f, 89.0f);
        PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
    }

    // Pan
    if (bPanMouseDown && !bIsFirstRightFrame && MouseDelta.SizeSquared() > 0.01f)
    {
        FVector RightVector = Camera->GetRightVector();
        FVector UpVector = Camera->GetUpVector();
        FVector NewLocation = GetActorLocation()
            - RightVector * MouseDelta.X * PanSpeed * 0.5f
            + UpVector * MouseDelta.Y * PanSpeed * 0.5f;
        SetActorLocation(NewLocation);
    }

    // Clear "first frame" flags
    bIsFirstLeftFrame = false;
    bIsFirstRightFrame = false;

    // Update tracking
    LastMousePosition = CurrentMousePos;
    bLeftMouseWasDown = bRotateMouseDown;
    bRightMouseWasDown = bPanMouseDown;

    // Zoom
    float WheelAxis = PC->GetInputAnalogKeyState(EKeys::MouseWheelAxis);
    if (FMath::Abs(WheelAxis) > 0.01f)
    {
        OrbitDistance = FMath::Clamp(OrbitDistance - WheelAxis * ZoomSpeed, MinZoomDistance, MaxZoomDistance);
        Camera->SetRelativeLocation(FVector(-OrbitDistance, 0, 0));
    }

    // Fit to screen with spacebar
    static bool bSpaceWasPressed = false;
    if (PC->IsInputKeyDown(EKeys::SpaceBar))
    {
        if (!bSpaceWasPressed)
        {
            FitActorToScreen();
            bSpaceWasPressed = true;
        }
    }
    else
    {
        bSpaceWasPressed = false;
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

void APDBCameraComponent::FitActorToScreen()
{
    // Find the PDBViewer actor
    TArray<AActor*> AllActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), AActor::StaticClass(), AllActors);

    UE_LOG(LogTemp, Warning, TEXT("Looking for PDBViewer actor among %d actors"), AllActors.Num());

    AActor* PDBViewerActor = nullptr;
    for (AActor* Actor : AllActors)
    {
        if (Actor == this)
            continue; // 🚫 Skip the camera itself

        FString ActorName = Actor->GetName();
        UE_LOG(LogTemp, Log, TEXT("Found actor: %s"), *ActorName);

        // Find actor with "PDB" in the name
        if (ActorName.Contains(TEXT("PDB"), ESearchCase::IgnoreCase))
        {
            PDBViewerActor = Actor;
            UE_LOG(LogTemp, Warning, TEXT("Found PDBViewer: %s"), *ActorName);
            break;
        }
    }

    // Fallback: find actor with most primitive components (likely the molecule)
    if (!PDBViewerActor)
    {
        UE_LOG(LogTemp, Error, TEXT("Could not find PDBViewer actor by name; searching by component count..."));
        int32 MaxComponents = 0;
        for (AActor* Actor : AllActors)
        {
            if (Actor == this || Actor->IsA(APawn::StaticClass()))
                continue;

            TArray<UPrimitiveComponent*> Components;
            Actor->GetComponents<UPrimitiveComponent>(Components);
            if (Components.Num() > MaxComponents)
            {
                MaxComponents = Components.Num();
                PDBViewerActor = Actor;
            }
        }

        if (PDBViewerActor)
        {
            UE_LOG(LogTemp, Warning, TEXT("Using actor with most components (%d): %s"),
                MaxComponents, *PDBViewerActor->GetName());
        }
    }

    if (!PDBViewerActor)
    {
        UE_LOG(LogTemp, Error, TEXT("Still could not find any suitable actor to focus on!"));
        return;
    }

    // Compute bounds
    FBox BoundingBox(ForceInit);
    bool bFoundAnyBounds = false;

    TArray<UPrimitiveComponent*> PrimitiveComponents;
    PDBViewerActor->GetComponents<UPrimitiveComponent>(PrimitiveComponents);

    UE_LOG(LogTemp, Warning, TEXT("PDBViewer has %d primitive components"), PrimitiveComponents.Num());

    for (UPrimitiveComponent* Primitive : PrimitiveComponents)
    {
        if (Primitive && Primitive->IsVisible())
        {
            FBox ComponentBounds = Primitive->Bounds.GetBox();
            if (ComponentBounds.IsValid)
            {
                BoundingBox += ComponentBounds;
                bFoundAnyBounds = true;
            }
        }
    }

    if (!bFoundAnyBounds)
    {
        // Fallback to actor location only
        SetActorLocation(PDBViewerActor->GetActorLocation());
        UE_LOG(LogTemp, Log, TEXT("Centered on PDBViewer actor location: %s"),
            *PDBViewerActor->GetActorLocation().ToString());

        CurrentYaw = 0.0f;
        CurrentPitch = -30.0f;
        PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
        return;
    }

    // Center and reset orientation
    FVector Center = BoundingBox.GetCenter();

    CurrentYaw = 0.0f;
    CurrentPitch = -30.0f;
    PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));

    SetActorLocation(Center);

    UE_LOG(LogTemp, Log, TEXT("Centered on PDBViewer at %s, keeping zoom distance %.1f"),
        *Center.ToString(), OrbitDistance);
}

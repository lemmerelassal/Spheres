#include "PDBCameraComponent.h"
#include "Camera/CameraComponent.h"
#include "Components/SceneComponent.h"
#include "GameFramework/PlayerController.h"
#include "Components/InputComponent.h"
#include "Framework/Application/SlateApplication.h"
#include "Kismet/GameplayStatics.h"
#include "PDBViewer.h"
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
    MinZoomDistance = 10.0f;
    MaxZoomDistance = 50000000.0f; // Much larger max distance

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

    // Fit to screen with spacebar - just center and reset rotation, keep current zoom
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
    // First try to center on the currently visible ligand (if any)
    TArray<AActor*> FoundViewers;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundViewers);

    for (AActor* Actor : FoundViewers)
    {
        APDBViewer* Viewer = Cast<APDBViewer>(Actor);
        if (!Viewer) continue;

        FLigandInfo* BestInfo = Viewer->GetVisibleLigandInfo();

        if (BestInfo && (BestInfo->AtomMeshes.Num() + BestInfo->BondMeshes.Num()) > 0)
        {
            // Compute bounding box from visible components of this ligand
            FBox BoundingBox(ForceInit);
            bool bFoundAnyBounds = false;

            for (UStaticMeshComponent* Mesh : BestInfo->AtomMeshes)
            {
                if (Mesh && Mesh->IsVisible())
                {
                    FBox CompBounds = Mesh->Bounds.GetBox();
                    if (CompBounds.IsValid)
                    {
                        BoundingBox += CompBounds;
                        bFoundAnyBounds = true;
                    }
                }
            }

            for (UStaticMeshComponent* Mesh : BestInfo->BondMeshes)
            {
                if (Mesh && Mesh->IsVisible())
                {
                    FBox CompBounds = Mesh->Bounds.GetBox();
                    if (CompBounds.IsValid)
                    {
                        BoundingBox += CompBounds;
                        bFoundAnyBounds = true;
                    }
                }
            }

            if (bFoundAnyBounds)
            {
                FVector Center = BoundingBox.GetCenter();
                SetActorLocation(Center);
                UE_LOG(LogTemp, Log, TEXT("Centered on visible ligand at %s, zoom=%.1f"),
                    *Center.ToString(), OrbitDistance);

                // Reset rotation to default view
                CurrentYaw = 0.0f;
                CurrentPitch = -30.0f;
                PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
                return;
            }
            else if (BestInfo->AtomMeshes.Num() > 0 && BestInfo->AtomMeshes[0])
            {
                // Fallback to a component location if no valid bounds
                SetActorLocation(BestInfo->AtomMeshes[0]->GetComponentLocation());
                CurrentYaw = 0.0f;
                CurrentPitch = -30.0f;
                PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
                return;
            }
        }
    }

    // Fallback: find actor with the most primitive components (that's our molecule)
    TArray<AActor*> AllActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), AActor::StaticClass(), AllActors);
    
    AActor* MoleculeActor = nullptr;
    int32 MaxComponents = 0;
    
    for (AActor* Actor : AllActors)
    {
        // Skip self and pawns
        if (Actor == this || Actor->IsA(APawn::StaticClass()))
            continue;
            
        TArray<UStaticMeshComponent*> MeshComponents;
        Actor->GetComponents<UStaticMeshComponent>(MeshComponents);
        
        // The molecule will have tons of sphere mesh components
        if (MeshComponents.Num() > MaxComponents && MeshComponents.Num() > 100)
        {
            MaxComponents = MeshComponents.Num();
            MoleculeActor = Actor;
        }
    }
    
    if (!MoleculeActor)
    {
        UE_LOG(LogTemp, Error, TEXT("Could not find molecule actor!"));
        return;
    }
    
    UE_LOG(LogTemp, Warning, TEXT("Found molecule actor '%s' with %d components"), 
        *MoleculeActor->GetName(), MaxComponents);
    
    // Get bounds of all the mesh components
    FBox BoundingBox(ForceInit);
    bool bFoundAnyBounds = false;
    
    TArray<UStaticMeshComponent*> MeshComponents;
    MoleculeActor->GetComponents<UStaticMeshComponent>(MeshComponents);
    
    for (UStaticMeshComponent* Mesh : MeshComponents)
    {
        if (Mesh && Mesh->IsVisible())
        {
            FBox ComponentBounds = Mesh->Bounds.GetBox();
            if (ComponentBounds.IsValid)
            {
                BoundingBox += ComponentBounds;
                bFoundAnyBounds = true;
            }
        }
    }
    
    if (!bFoundAnyBounds)
    {
        // Fallback: just use the actor's location
        SetActorLocation(MoleculeActor->GetActorLocation());
        UE_LOG(LogTemp, Warning, TEXT("No bounds found, centered on actor location: %s"), 
            *MoleculeActor->GetActorLocation().ToString());
    }
    else
    {
        // Get center
        FVector Center = BoundingBox.GetCenter();
        SetActorLocation(Center);
        UE_LOG(LogTemp, Log, TEXT("Centered on molecule at %s, zoom=%.1f"), 
            *Center.ToString(), OrbitDistance);
    }
    
    // Reset rotation to default view
    CurrentYaw = 0.0f;
    CurrentPitch = -30.0f;
    PivotPoint->SetRelativeRotation(FRotator(CurrentPitch, CurrentYaw, 0));
}
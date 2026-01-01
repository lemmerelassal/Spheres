#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Pawn.h"
#include "PDBCameraComponent.generated.h"

UCLASS()
class SPHERES_API APDBCameraComponent : public APawn
{
    GENERATED_BODY()

public:
    APDBCameraComponent();

protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;
    virtual void SetupPlayerInputComponent(class UInputComponent* PlayerInputComponent) override;

    // Set the target actor to orbit around
    UFUNCTION(BlueprintCallable, Category = "Camera")
    void SetTargetActor(AActor* Target);

    // Fit the target actor to screen
    void FitActorToScreen();

private:
    // Camera components
    UPROPERTY(VisibleAnywhere)
    class USceneComponent* PivotPoint;

    UPROPERTY(VisibleAnywhere)
    class UCameraComponent* Camera;

    // Input tracking
    FVector2D LastMousePosition;
    FVector2D MousePositionOnPress;
    bool bLeftMouseWasDown;
    bool bRightMouseWasDown;
    bool bIsFirstLeftFrame;
    bool bIsFirstRightFrame;

    // Camera parameters
    UPROPERTY(EditAnywhere, Category = "Camera")
    float OrbitDistance;

    UPROPERTY(EditAnywhere, Category = "Camera")
    float RotationSpeed;

    UPROPERTY(EditAnywhere, Category = "Camera")
    float PanSpeed;

    UPROPERTY(EditAnywhere, Category = "Camera")
    float ZoomSpeed;

    UPROPERTY(EditAnywhere, Category = "Camera")
    float MinZoomDistance;

    UPROPERTY(EditAnywhere, Category = "Camera")
    float MaxZoomDistance;

    // Current rotation
    float CurrentYaw;
    float CurrentPitch;
};
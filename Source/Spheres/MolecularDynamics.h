// MolecularDynamics.h - Molecular Dynamics Engine for PDB Viewer
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "MolecularDynamics.generated.h"

// Forward declarations
class APDBViewer;
class UStaticMeshComponent;

USTRUCT(BlueprintType)
struct FAtomState
{
    GENERATED_BODY()
    
    UPROPERTY()
    FVector Position;
    
    UPROPERTY()
    FVector Velocity;
    
    UPROPERTY()
    FVector Force;
    
    UPROPERTY()
    float Mass;
    
    UPROPERTY()
    FString Element;
    
    UPROPERTY()
    UStaticMeshComponent* MeshComponent;
    
    FAtomState()
        : Position(FVector::ZeroVector)
        , Velocity(FVector::ZeroVector)
        , Force(FVector::ZeroVector)
        , Mass(1.0f)
        , Element(TEXT("C"))
        , MeshComponent(nullptr)
    {}
};

USTRUCT(BlueprintType)
struct FBondConstraint
{
    GENERATED_BODY()
    
    UPROPERTY()
    int32 Atom1Index;
    
    UPROPERTY()
    int32 Atom2Index;
    
    UPROPERTY()
    float RestLength;
    
    UPROPERTY()
    float SpringConstant;
    
    UPROPERTY()
    int32 BondOrder;
    
    FBondConstraint()
        : Atom1Index(-1)
        , Atom2Index(-1)
        , RestLength(100.0f)
        , SpringConstant(500.0f)
        , BondOrder(1)
    {}
};

UENUM(BlueprintType)
enum class EMDIntegrator : uint8
{
    Euler UMETA(DisplayName = "Euler"),
    Verlet UMETA(DisplayName = "Verlet"),
    LeapFrog UMETA(DisplayName = "Leap-Frog")
};

UCLASS()
class SPHERES_API AMolecularDynamics : public AActor
{
    GENERATED_BODY()

public:
    AMolecularDynamics();
    
    virtual void Tick(float DeltaTime) override;
    
    // Simulation control
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void StartSimulation();
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void StopSimulation();
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void ResetSimulation();
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void StepSimulation(float DeltaTime);
    
    // Initialize from PDB Viewer
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void InitializeFromViewer(APDBViewer* Viewer);
    
    // Parameter setters
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetTimeStep(float NewTimeStep) { TimeStep = NewTimeStep; }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetTemperature(float NewTemperature) { TargetTemperature = NewTemperature; }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetDamping(float NewDamping) { DampingFactor = FMath::Clamp(NewDamping, 0.0f, 1.0f); }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetIntegrator(EMDIntegrator NewIntegrator) { Integrator = NewIntegrator; }
    
    // Getters
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    bool IsSimulating() const { return bIsSimulating; }
    
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    float GetCurrentEnergy() const { return TotalEnergy; }
    
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    float GetKineticEnergy() const { return KineticEnergy; }
    
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    float GetPotentialEnergy() const { return PotentialEnergy; }
    
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    int32 GetAtomCount() const { return Atoms.Num(); }

protected:
    virtual void BeginPlay() override;
    
    // Simulation parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float TimeStep = 0.001f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float TargetTemperature = 300.0f; // Kelvin
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float DampingFactor = 0.95f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    EMDIntegrator Integrator = EMDIntegrator::Verlet;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    bool bUseLennardJones = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    bool bUseElectrostatics = false;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float LJEpsilon = 0.1f; // Lennard-Jones well depth
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float LJSigma = 350.0f; // Lennard-Jones distance parameter
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    float CutoffDistance = 1000.0f; // Interaction cutoff
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    bool bConstrainBonds = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Dynamics")
    int32 ConstraintIterations = 5;

private:
    // Simulation state
    UPROPERTY()
    TArray<FAtomState> Atoms;
    
    UPROPERTY()
    TArray<FBondConstraint> Bonds;
    
    UPROPERTY()
    TArray<FVector> InitialPositions;
    
    UPROPERTY()
    APDBViewer* ViewerReference;
    
    bool bIsSimulating = false;
    float SimulationTime = 0.0f;
    float TotalEnergy = 0.0f;
    float KineticEnergy = 0.0f;
    float PotentialEnergy = 0.0f;
    
    // Previous positions for Verlet integration
    TArray<FVector> PreviousPositions;
    
    // Integration methods
    void IntegrateEuler(float DeltaTime);
    void IntegrateVerlet(float DeltaTime);
    void IntegrateLeapFrog(float DeltaTime);
    
    // Force calculations
    void CalculateForces();
    void CalculateBondForces();
    void CalculateLennardJonesForces();
    void CalculateElectrostaticForces();
    
    // Constraint enforcement
    void EnforceConstraints();
    
    // Energy calculations
    void UpdateEnergies();
    
    // Temperature control
    void ApplyThermostat();
    
    // Helper functions
    float GetAtomMass(const FString& Element) const;
    void UpdateMeshPositions();
};
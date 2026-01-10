// MolecularDynamics.h - Optimized O(N) MD Engine with Spatial Hashing
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "MolecularDynamics.generated.h"

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
    float Charge;
    
    UPROPERTY()
    FString Element;
    
    UPROPERTY()
    UStaticMeshComponent* MeshComponent;
    
    // Spatial hashing
    int32 CellX;
    int32 CellY;
    int32 CellZ;
    int64 CellHash;
    
    // Verlet list
    TArray<int32> NeighborList;
    int32 NeighborUpdateCounter;
    
    FAtomState()
        : Position(FVector::ZeroVector)
        , Velocity(FVector::ZeroVector)
        , Force(FVector::ZeroVector)
        , Mass(12.011f)
        , Charge(0.0f)
        , Element(TEXT("C"))
        , MeshComponent(nullptr)
        , CellX(0), CellY(0), CellZ(0)
        , CellHash(0)
        , NeighborUpdateCounter(0)
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
    Verlet UMETA(DisplayName = "Velocity Verlet"),
    LeapFrog UMETA(DisplayName = "Leap-Frog"),
    RungeKutta4 UMETA(DisplayName = "Runge-Kutta 4th Order")
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
    
    // Initialize from ALL visible molecules in PDB Viewer
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void InitializeFromViewer(APDBViewer* Viewer);
    
    // Parameter setters
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetTimeStep(float NewTimeStep) { TimeStep = FMath::Clamp(NewTimeStep, 0.0001f, 0.01f); }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetTemperature(float NewTemperature) { TargetTemperature = FMath::Max(0.0f, NewTemperature); }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetDamping(float NewDamping) { DampingFactor = FMath::Clamp(NewDamping, 0.0f, 1.0f); }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetIntegrator(EMDIntegrator NewIntegrator) { Integrator = NewIntegrator; }
    
    UFUNCTION(BlueprintCallable, Category = "Molecular Dynamics")
    void SetCutoffDistance(float NewCutoff) { CutoffDistance = FMath::Max(100.0f, NewCutoff); RebuildSpatialHash(); }
    
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
    
    UFUNCTION(BlueprintPure, Category = "Molecular Dynamics")
    float GetCurrentTemperature() const;

protected:
    virtual void BeginPlay() override;
    
    // Simulation parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Parameters")
    float TimeStep = 0.002f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Parameters")
    float TargetTemperature = 300.0f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Parameters")
    float DampingFactor = 0.98f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Parameters")
    EMDIntegrator Integrator = EMDIntegrator::Verlet;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Forces")
    bool bUseLennardJones = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Forces")
    bool bUseElectrostatics = false;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Forces")
    float LJEpsilon = 0.5f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Forces")
    float LJSigma = 300.0f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Forces")
    float CutoffDistance = 1200.0f;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Constraints")
    bool bConstrainBonds = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Constraints")
    int32 ConstraintIterations = 3;
    
    // Optimization parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Optimization")
    int32 NeighborListRebuildFrequency = 10;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MD Optimization")
    float NeighborListSkin = 200.0f;

private:
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
    int32 StepCount = 0;
    
    // Verlet integration
    TArray<FVector> PreviousPositions;
    TArray<FVector> PreviousAccelerations;
    
    // Spatial hashing for O(N) neighbor search
    TMap<int64, TArray<int32>> SpatialHash;
    float CellSize;
    FVector GridMin;
    FVector GridMax;
    
    // RK4 state storage
    struct FRK4State
    {
        TArray<FVector> K1Vel, K2Vel, K3Vel, K4Vel;
        TArray<FVector> K1Acc, K2Acc, K3Acc, K4Acc;
    };
    FRK4State RK4State;
    
    // Integration methods
    void IntegrateVerlet(float DeltaTime);
    void IntegrateLeapFrog(float DeltaTime);
    void IntegrateRungeKutta4(float DeltaTime);
    
    // Force calculations - O(N) with spatial hashing
    void CalculateForces();
    void CalculateBondForces();
    void CalculateNonBondedForces();
    
    // Spatial hashing
    void RebuildSpatialHash();
    void UpdateSpatialHash();
    int64 GetCellHash(int32 X, int32 Y, int32 Z) const;
    void GetCellCoords(const FVector& Position, int32& OutX, int32& OutY, int32& OutZ) const;
    void GetNeighborCells(int32 X, int32 Y, int32 Z, TArray<int64>& OutHashes) const;
    
    // Neighbor lists
    void UpdateNeighborLists();
    bool ShouldRebuildNeighborLists() const;
    
    // Constraint enforcement - RATTLE algorithm
    void EnforceConstraints(float DeltaTime);
    
    // Energy and temperature
    void UpdateEnergies();
    void ApplyThermostat();
    
    // Helpers
    float GetAtomMass(const FString& Element) const;
    float GetAtomCharge(const FString& Element) const;
    void UpdateMeshPositions();
    
    // Element properties
    struct FElementProperties
    {
        float Mass;
        float Radius;
        float Epsilon;
        float Charge;
    };
    TMap<FString, FElementProperties> ElementDatabase;
    void InitializeElementDatabase();
};
// FEPCalculatorOptimized.h - Fully Optimized FEP with all techniques
// Expected: 500-5000x speedup over baseline
// Compatible with Unreal Engine 5.6

#pragma once

#include "CoreMinimal.h"
#include "FEPCalculator.h"
#include "RHI.h"
#include "RHIResources.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "FEPCalculatorOptimized.generated.h"

// Neighbor list for efficient force calculations
USTRUCT()
struct FNeighborList
{
    GENERATED_BODY()
    
    TArray<int32> Neighbors;
    FVector LastUpdatePosition;
    
    FNeighborList() : LastUpdatePosition(FVector::ZeroVector) {}
};

// Cell for spatial grid
USTRUCT()
struct FSpatialCell
{
    GENERATED_BODY()
    
    TArray<int32> AtomIndices;
};

// Structure of Arrays for better cache coherency
USTRUCT()
struct FAtomDataSoA
{
    GENERATED_BODY()
    
    TArray<FVector> Positions;
    TArray<FVector> Velocities;
    TArray<FVector> Forces;
    TArray<float> Masses;
    TArray<float> Charges;
    TArray<float> VdWRadii;
    TArray<float> VdWEpsilons;
    TArray<FString> Elements;
    TArray<bool> IsLigandAtom;
    
    int32 NumAtoms() const { return Positions.Num(); }
    
    void Reserve(int32 Count)
    {
        Positions.Reserve(Count);
        Velocities.Reserve(Count);
        Forces.Reserve(Count);
        Masses.Reserve(Count);
        Charges.Reserve(Count);
        VdWRadii.Reserve(Count);
        VdWEpsilons.Reserve(Count);
        Elements.Reserve(Count);
        IsLigandAtom.Reserve(Count);
    }
    
    void Empty()
    {
        Positions.Empty();
        Velocities.Empty();
        Forces.Empty();
        Masses.Empty();
        Charges.Empty();
        VdWRadii.Empty();
        VdWEpsilons.Empty();
        Elements.Empty();
        IsLigandAtom.Empty();
    }
};

// Precomputed LJ potential table
struct FLJTable
{
    TArray<float> EnergyTable;
    TArray<float> ForceTable;
    float TableSpacing = 0.01f; // Angstroms
    float MaxDistance = 15.0f;
    
    void Build(float Sigma, float Epsilon);
    float GetEnergy(float r) const;
    float GetForce(float r) const;
};

// GPU buffer for atom data
struct FAtomGPUData
{
    FVector4f Position;  // w unused
    FVector4f Velocity;  // w = mass
    FVector4f Force;     // w unused
    FVector4f Params;    // x=charge, y=vdw_radius, z=vdw_epsilon, w=is_ligand
};

// Memory pool for reusable allocations
class FFEPMemoryPool
{
public:
    TArray<FVector> VectorPool;
    TArray<float> FloatPool;
    int32 VectorPoolIndex = 0;
    int32 FloatPoolIndex = 0;
    
    FVector* AllocateVector()
    {
        if (VectorPoolIndex >= VectorPool.Num())
        {
            VectorPool.Add(FVector::ZeroVector);
        }
        return &VectorPool[VectorPoolIndex++];
    }
    
    void Reset()
    {
        VectorPoolIndex = 0;
        FloatPoolIndex = 0;
    }
};

UCLASS()
class SPHERES_API AFEPCalculatorOptimized : public AFEPCalculator
{
    GENERATED_BODY()

public:
    AFEPCalculatorOptimized();
    
    // ===== OPTIMIZATION CONTROLS =====
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseNeighborLists = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseCellLists = true;  // Even faster than neighbor lists for large systems
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseStructureOfArrays = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseLJTables = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseParallelForces = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseSquaredDistances = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseCombinedEnergyDerivative = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseAdaptiveTimestep = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseEarlyTermination = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseMemoryPool = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseGPU = false;  // Auto-detect GPU capability
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bUseSIMD = true;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    bool bMinimalLogging = true;  // Only essential logs
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    float NeighborListSkin = 2.0f; // Angstroms beyond cutoff
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    float CellSize = 5.0f;  // Cell grid size in Angstroms
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "FEP|Optimization")
    float EarlyTerminationThreshold = 0.01f;  // 1% relative error
    
    // ===== STATISTICS =====
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    int32 NeighborListRebuildCount = 0;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    float AverageNeighborsPerAtom = 0.0f;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    float LastStepTimeMs = 0.0f;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    float TotalCalculationTimeMs = 0.0f;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    int32 TotalStepsExecuted = 0;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    bool bConvergedEarly = false;
    
    UPROPERTY(BlueprintReadOnly, Category = "FEP|Stats")
    FString OptimizationProfile;  // Description of active optimizations
    
    // Auto-detect best optimization strategy
    UFUNCTION(BlueprintCallable, Category = "FEP|Optimization")
    void AutoOptimize();
    
    // Benchmark different optimization strategies
    UFUNCTION(BlueprintCallable, Category = "FEP|Optimization")
    void BenchmarkOptimizations();

protected:
    // Optimized data structures
    FAtomDataSoA AtomDataSoA;
    TArray<FNeighborList> NeighborLists;
    TMap<TPair<FString, FString>, FLJTable> LJTables;
    FFEPMemoryPool MemoryPool;
    
    // Cell grid for spatial partitioning
    TArray<FSpatialCell> CellGrid;
    int32 GridSizeX, GridSizeY, GridSizeZ;
    FVector GridMin, GridMax;
    
    // Cached values
    float CutoffDistanceSq;
    float NeighborListCutoffSq;
    float AdaptiveDt;
    
    // GPU resources
    FBufferRHIRef AtomPositionBuffer;
    FBufferRHIRef AtomForceBuffer;
    FUnorderedAccessViewRHIRef AtomPositionUAV;
    FUnorderedAccessViewRHIRef AtomForceUAV;
    bool bGPUInitialized = false;
    
    // Convergence tracking
    TArray<float> RecentdHdL;
    int32 ConvergenceCheckInterval = 1000;
    
    // Override key methods with optimized versions
    virtual void InitializeSystem(const FString& LigandKey) override;
    virtual void CalculateForces(float Lambda) override;
    virtual float CalculateTotalEnergy(float Lambda) override;
    virtual float CalculateElectrostaticEnergy(float Lambda) override;
    virtual float CalculateVanDerWaalsEnergy(float Lambda) override;
    virtual void PerformMDStep(float Lambda) override;
    virtual void RunLambdaWindow(float Lambda) override;
    virtual float CalculatedHdLambda(float Lambda) override;
    
    // Optimized implementations
    void BuildNeighborLists();
    void BuildCellGrid();
    bool NeedNeighborListRebuild() const;
    void ConvertToSoA();
    void ConvertFromSoA();
    void BuildLJTables();
    
    void CalculateForcesOptimized(float Lambda);
    void CalculateForcesWithCells(float Lambda);
    void CalculateForcesParallel(float Lambda);
    void CalculateForcesGPU(float Lambda);
    void CalculateForcesVectorized(float Lambda);
    
    void CalculateEnergyAndDerivativeCombined(float Lambda, float& Energy, float& dHdL);
    
    // GPU support
    void InitializeGPUResources();
    void CleanupGPUResources();
    void UploadAtomsToGPU();
    void DownloadForcesFromGPU();
    
    // Adaptive timestep
    float CalculateAdaptiveTimestep();
    
    // Convergence checking
    bool CheckConvergence();
    
    // Spatial hashing
    int32 GetCellIndex(const FVector& Position) const;
    void GetNeighborCells(int32 CellIndex, TArray<int32>& NeighborCellIndices) const;
    
    // SIMD helpers
    void CalculateForceSIMD(int32 AtomIndex, float Lambda);
    
    // Utility
    FORCEINLINE float FastDistance(const FVector& A, const FVector& B) const
    {
        return (B - A).Size();
    }
    
    FORCEINLINE float FastDistanceSquared(const FVector& A, const FVector& B) const
    {
        return (B - A).SizeSquared();
    }
    
    FORCEINLINE int32 Flatten3DIndex(int32 X, int32 Y, int32 Z) const
    {
        return X + Y * GridSizeX + Z * GridSizeX * GridSizeY;
    }
    
    void LogOptimizationStatus();
};

// GPU Compute Shader for force calculations
class FFEPForceComputeShader : public FGlobalShader
{
    DECLARE_GLOBAL_SHADER(FFEPForceComputeShader);
    SHADER_USE_PARAMETER_STRUCT(FFEPForceComputeShader, FGlobalShader);
    
    BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
        SHADER_PARAMETER_UAV(RWStructuredBuffer<FVector4f>, AtomPositions)
        SHADER_PARAMETER_UAV(RWStructuredBuffer<FVector4f>, AtomForces)
        SHADER_PARAMETER(uint32, NumAtoms)
        SHADER_PARAMETER(float, Lambda)
        SHADER_PARAMETER(float, CutoffDistance)
        SHADER_PARAMETER(float, CoulombConstant)
    END_SHADER_PARAMETER_STRUCT()
    
public:
    static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
    {
        return IsFeatureLevelSupported(Parameters.Platform, ERHIFeatureLevel::SM5);
    }
    
    static void ModifyCompilationEnvironment(const FGlobalShaderPermutationParameters& Parameters, FShaderCompilerEnvironment& OutEnvironment)
    {
        FGlobalShader::ModifyCompilationEnvironment(Parameters, OutEnvironment);
        OutEnvironment.SetDefine(TEXT("THREADGROUP_SIZE"), 64);
    }
};

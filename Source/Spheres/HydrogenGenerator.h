// HydrogenGenerator.h - Adds explicit hydrogens to molecular structures

#pragma once

#include "CoreMinimal.h"

struct FAtomBondInfo
{
    FString Element;
    FVector Position;
    TArray<int32> BondedAtomIndices;
    TArray<int32> BondOrders; // Track bond order for each bonded atom
    int32 BondCount = 0;
    int32 TotalBondOrder = 0; // Sum of all bond orders (single=1, double=2, triple=3)
};

class SPHERES_API FHydrogenGenerator
{
public:
    // Calculate and add hydrogens to a set of atoms
    // BondOrders array contains bond order for each bond pair (aligned with BondPairs)
    // Returns array of hydrogen positions and their bonded heavy atom index
    static TArray<TPair<FVector, int32>> GenerateHydrogens(
        const TArray<FVector>& AtomPositions,
        const TArray<FString>& AtomElements,
        const TArray<TPair<int32, int32>>& BondPairs,
        const TArray<int32>& BondOrders);

private:
    // Get ideal number of hydrogens based on valence and TOTAL bond order
    static int32 GetIdealHydrogenCount(const FString& Element, int32 TotalBondOrder);
    
    // Get ideal bond length for element-hydrogen bond
    static float GetIdealBondLength(const FString& Element);
    
    // Generate hydrogen positions for different geometries
    static TArray<FVector> GenerateHydrogensForAtom(
        const FVector& AtomPos,
        const FString& Element,
        const TArray<FVector>& BondedPositions,
        int32 NumHydrogens,
        int32 TotalBondOrder);
    
    // Geometry-specific functions
    static TArray<FVector> GenerateTetrahedralHydrogens(
        const FVector& Center,
        const TArray<FVector>& Existing,
        int32 NumH,
        float BondLength);
    
    static TArray<FVector> GenerateTrigonalHydrogens(
        const FVector& Center,
        const TArray<FVector>& Existing,
        int32 NumH,
        float BondLength);
    
    static TArray<FVector> GenerateLinearHydrogens(
        const FVector& Center,
        const TArray<FVector>& Existing,
        int32 NumH,
        float BondLength);
    
    // Helper to get perpendicular vector
    static FVector GetPerpendicularVector(const FVector& V);
};
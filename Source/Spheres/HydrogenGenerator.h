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
    /**
     * Generate hydrogen positions for a molecule
     * @param AtomPositions - Positions of all heavy atoms
     * @param AtomElements - Element symbols for each atom
     * @param BondPairs - Connectivity information (atom index pairs)
     * @param BondOrders - Bond order for each bond (1=single, 2=double, 3=triple)
     * @return Array of hydrogen positions paired with their parent atom index
     */
    static TArray<TPair<FVector, int32>> GenerateHydrogens(
        const TArray<FVector>& AtomPositions,
        const TArray<FString>& AtomElements,
        const TArray<TPair<int32, int32>>& BondPairs,
        const TArray<int32>& BondOrders
    );

private:
    /**
     * Get the typical valence (bonding capacity) for an element
     * @param Element - Element symbol (e.g., "C", "N", "O")
     * @return Typical valence, or 0 if unknown
     */
    static int32 GetTypicalValence(const FString& Element);


    /**
     * Generate 3D positions for hydrogens based on parent atom geometry
     * @param ParentPos - Position of the parent heavy atom
     * @param ParentIdx - Index of the parent atom
     * @param NumHydrogens - Number of hydrogens to generate
     * @param AllAtomPositions - All atom positions (for geometric context)
     * @param BondPairs - Bond connectivity information
     * @param ParentElement - Element symbol of parent atom
     * @return Array of hydrogen positions
     */
    static TArray<FVector> GenerateHydrogenPositions(
        const FVector& ParentPos,
        int32 ParentIdx,
        int32 NumHydrogens,
        const TArray<FVector>& AllAtomPositions,
        const TArray<TPair<int32, int32>>& BondPairs,
        const FString& ParentElement
    );

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
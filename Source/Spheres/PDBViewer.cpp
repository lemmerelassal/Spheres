#include "PDBViewer.h"
#include "HttpModule.h"
#include "Interfaces/IHttpResponse.h"
#include "DrawDebugHelpers.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Engine/StaticMesh.h"
#include "Engine/World.h"
#include "Components/StaticMeshComponent.h"
#include "Components/SceneComponent.h"

// ----------------------------
// Constructor
// ----------------------------
APDBViewer::APDBViewer()
{
    PrimaryActorTick.bCanEverTick = false;

    // Root component
    USceneComponent* RootComp = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = RootComp;

    // Load default sphere mesh
    static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereMeshObj(TEXT("/Engine/BasicShapes/Sphere.Sphere"));
    if (SphereMeshObj.Succeeded()) SphereMeshAsset = SphereMeshObj.Object;

    // Load default cylinder mesh (for bonds, optional)
    static ConstructorHelpers::FObjectFinder<UStaticMesh> CylinderMeshObj(TEXT("/Engine/BasicShapes/Cylinder.Cylinder"));
    if (CylinderMeshObj.Succeeded()) CylinderMeshAsset = CylinderMeshObj.Object;

    // Load default material
    static ConstructorHelpers::FObjectFinder<UMaterialInterface> SphereMatObj(TEXT("/Engine/BasicShapes/BasicShapeMaterial.BasicShapeMaterial"));
    if (SphereMatObj.Succeeded()) SphereMaterialAsset = SphereMatObj.Object;
}

// ----------------------------
// BeginPlay
// ----------------------------
void APDBViewer::BeginPlay()
{
    Super::BeginPlay();

    // Example: fetch a PDB structure
    FetchAndDisplayStructure(TEXT("5ENB"));
}

// ----------------------------
// Fetch PDB or CIF
// ----------------------------
void APDBViewer::FetchAndDisplayStructure(const FString& PDB_ID)
{
    FString PDB_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.pdb"), *PDB_ID);
    FString CIF_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.cif"), *PDB_ID);

    UE_LOG(LogTemp, Log, TEXT("Fetching PDB file: %s"), *PDB_URL);

    FetchFileAsync(PDB_URL, [this, CIF_URL](bool bSuccess, const FString& Content)
    {
        if (bSuccess)
        {
            ParsePDB(Content);
        }
        else
        {
            UE_LOG(LogTemp, Warning, TEXT("PDB not found. Trying mmCIF..."));
            FetchFileAsync(CIF_URL, [this](bool bSuccess2, const FString& Content2)
            {
                if (bSuccess2)
                {
                    ParseMMCIF(Content2);
                }
                else
                {
                    UE_LOG(LogTemp, Error, TEXT("Failed to fetch both PDB and mmCIF files."));
                }
            });
        }
    });
}

// ----------------------------
// Async HTTP fetch
// ----------------------------
void APDBViewer::FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback)
{
    TSharedRef<IHttpRequest, ESPMode::ThreadSafe> Request = FHttpModule::Get().CreateRequest();
    Request->SetURL(URL);
    Request->SetVerb(TEXT("GET"));

    Request->OnProcessRequestComplete().BindLambda([Callback](FHttpRequestPtr Req, FHttpResponsePtr Response, bool bWasSuccessful)
    {
        if (bWasSuccessful && Response.IsValid() && Response->GetResponseCode() == 200)
            Callback(true, Response->GetContentAsString());
        else
            Callback(false, TEXT(""));
    });

    Request->ProcessRequest();
}

// ----------------------------
// Parse PDB (ATOM only, centered + scaled)
// ----------------------------
void APDBViewer::ParsePDB(const FString& FileContent)
{
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    TArray<FVector> AtomPositions;
    TArray<FString> AtomElements;

    for (const FString& Line : Lines)
    {
        if (!Line.StartsWith("ATOM")) continue;

        // Coordinates
        float X = FCString::Atof(*Line.Mid(30, 8));
        float Y = FCString::Atof(*Line.Mid(38, 8));
        float Z = FCString::Atof(*Line.Mid(46, 8));

        AtomPositions.Add(FVector(X, Y, Z));

        // Element symbol
        FString Element = Line.Mid(76, 2).TrimStartAndEnd().ToUpper();
        if (Element.IsEmpty())
        {
            FString AtomName = Line.Mid(12, 4).TrimStartAndEnd();
            Element = AtomName.Left(1).ToUpper();
        }
        AtomElements.Add(Element);
    }

    if (AtomPositions.Num() == 0) return;

    // Compute centroid
    FVector Centroid = FVector::ZeroVector;
    for (const FVector& Pos : AtomPositions) Centroid += Pos;
    Centroid /= AtomPositions.Num();

    // Draw spheres (scale + center)
    const float ScaleFactor = 10.f; // 1 Å = 10 Unreal units
    TArray<UStaticMeshComponent*> AtomMeshes;

    for (int32 i = 0; i < AtomPositions.Num(); i++)
    {
        FVector Pos = (AtomPositions[i] - Centroid) * ScaleFactor;
        DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(AtomElements[i]), GetRootComponent(), AtomMeshes);
    }

    UE_LOG(LogTemp, Log, TEXT("Parsed %d ATOM entries from PDB."), AtomPositions.Num());
}

// ----------------------------
// Parse mmCIF (ATOM only, centered + scaled)
// ----------------------------
void APDBViewer::ParseMMCIF(const FString& FileContent)
{
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    TArray<FString> Headers;
    TArray<TArray<FString>> AtomTable;
    int XIdx = -1, YIdx = -1, ZIdx = -1, AtomIdx = -1, ElementIdx = -1, GroupIdx = -1;
    bool bInLoop = false;

    for (const FString& Line : Lines)
    {
        if (Line.StartsWith("loop_"))
        {
            bInLoop = true;
            Headers.Empty();
            continue;
        }

        if (bInLoop && Line.StartsWith("_atom_site."))
        {
            Headers.Add(Line);
            if (Line.Contains("Cartn_x")) XIdx = Headers.Num() - 1;
            if (Line.Contains("Cartn_y")) YIdx = Headers.Num() - 1;
            if (Line.Contains("Cartn_z")) ZIdx = Headers.Num() - 1;
            if (Line.Contains("label_atom_id")) AtomIdx = Headers.Num() - 1;
            if (Line.Contains("type_symbol")) ElementIdx = Headers.Num() - 1;
            if (Line.Contains("group_PDB")) GroupIdx = Headers.Num() - 1;
            continue;
        }

        if (bInLoop && !Line.StartsWith("_"))
        {
            TArray<FString> Tokens;
            Line.ParseIntoArrayWS(Tokens);
            if (Tokens.Num() > FMath::Max3(XIdx, YIdx, ZIdx))
                AtomTable.Add(Tokens);
        }
    }

    TArray<FVector> AtomPositions;
    TArray<FString> AtomElements;

    for (const auto& Row : AtomTable)
    {
        if (XIdx < 0 || YIdx < 0 || ZIdx < 0 || AtomIdx < 0) continue;
        if (GroupIdx >= 0 && Row[GroupIdx] != "ATOM") continue; // ATOM only

        float X = FCString::Atof(*Row[XIdx]);
        float Y = FCString::Atof(*Row[YIdx]);
        float Z = FCString::Atof(*Row[ZIdx]);
        AtomPositions.Add(FVector(X, Y, Z));

        FString Element = (ElementIdx >= 0) ? Row[ElementIdx].ToUpper() : Row[AtomIdx].Left(1).ToUpper();
        AtomElements.Add(Element);
    }

    if (AtomPositions.Num() == 0) return;

    // Compute centroid
    FVector Centroid = FVector::ZeroVector;
    for (const FVector& Pos : AtomPositions) Centroid += Pos;
    Centroid /= AtomPositions.Num();

    const float ScaleFactor = 10.f;
    TArray<UStaticMeshComponent*> AtomMeshes;

    for (int32 i = 0; i < AtomPositions.Num(); i++)
    {
        FVector Pos = (AtomPositions[i] - Centroid) * ScaleFactor;
        DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(AtomElements[i]), GetRootComponent(), AtomMeshes);
    }

    UE_LOG(LogTemp, Log, TEXT("Parsed %d ATOM entries from mmCIF."), AtomPositions.Num());
}

// ----------------------------
// Draw sphere
// ----------------------------
void APDBViewer::DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    if (!SphereMeshAsset || !SphereMaterialAsset || !Parent) return;

    UStaticMeshComponent* Sphere = NewObject<UStaticMeshComponent>(this);
    Sphere->RegisterComponent();
    Sphere->AttachToComponent(Parent, FAttachmentTransformRules::KeepRelativeTransform);
    Sphere->SetStaticMesh(SphereMeshAsset);
    Sphere->SetRelativeLocation(FVector(x, y, z));
    Sphere->SetWorldScale3D(FVector(0.1f));

    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue(FName("Color"), Color);
    Mat->SetScalarParameterValue(FName("EmissiveIntensity"), 5.f);
    Sphere->SetMaterial(0, Mat);

    OutArray.Add(Sphere);
}

// ----------------------------
// Element colors
// ----------------------------
FLinearColor APDBViewer::ElementColor(const FString& Element)
{
    if (Element.Equals("C", ESearchCase::IgnoreCase)) return FLinearColor(0.1f, 0.1f, 0.1f);
    if (Element.Equals("O", ESearchCase::IgnoreCase)) return FLinearColor::Red;
    if (Element.Equals("H", ESearchCase::IgnoreCase) || Element.Equals("D", ESearchCase::IgnoreCase)) return FLinearColor::White;
    if (Element.Equals("N", ESearchCase::IgnoreCase)) return FLinearColor::Blue;
    if (Element.Equals("S", ESearchCase::IgnoreCase)) return FLinearColor::Yellow;
    if (Element.Equals("CL", ESearchCase::IgnoreCase)) return FLinearColor(0.0f, 1.0f, 0.0f);
    if (Element.Equals("P", ESearchCase::IgnoreCase)) return FLinearColor(1.0f, 0.5f, 0.0f);
    if (Element.Equals("F", ESearchCase::IgnoreCase)) return FLinearColor(0.0f, 1.0f, 0.0f);
    if (Element.Equals("BR", ESearchCase::IgnoreCase)) return FLinearColor(0.6f, 0.2f, 0.2f);
    if (Element.Equals("I", ESearchCase::IgnoreCase)) return FLinearColor(0.4f, 0.0f, 0.8f);
    if (Element.Equals("FE", ESearchCase::IgnoreCase)) return FLinearColor(0.8f, 0.4f, 0.0f);
    if (Element.Equals("MG", ESearchCase::IgnoreCase)) return FLinearColor(0.0f, 0.8f, 0.0f);
    if (Element.Equals("ZN", ESearchCase::IgnoreCase)) return FLinearColor(0.5f, 0.5f, 0.5f);
    if (Element.Equals("CA", ESearchCase::IgnoreCase)) return FLinearColor(0.2f, 0.6f, 1.0f);
    if (Element.Equals("NA", ESearchCase::IgnoreCase)) return FLinearColor(0.0f, 0.0f, 1.0f);
    if (Element.Equals("K", ESearchCase::IgnoreCase)) return FLinearColor(0.5f, 0.0f, 1.0f);
    if (Element.Equals("CU", ESearchCase::IgnoreCase)) return FLinearColor(1.0f, 0.5f, 0.0f);
    if (Element.Equals("B", ESearchCase::IgnoreCase)) return FLinearColor(1.0f, 0.7f, 0.7f);
    return FLinearColor::Gray;
}

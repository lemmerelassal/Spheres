// MoleculeListWidget.cpp
#include "MoleculeListWidget.h"
#include "PDBViewer.h"
#include "Components/ListView.h"
#include "Kismet/GameplayStatics.h"

void UMoleculeListWidget::NativeConstruct()
{
    Super::NativeConstruct();

    // Automatically populate when widget is created
    AutoPopulate();
}

void UMoleculeListWidget::PopulateList(APDBViewer* Viewer)
{
    if (!Viewer || !MoleculeListView)
        return;

    CachedViewer = Viewer;

    // Unbind previous delegates
    if (CachedViewer)
    {
        CachedViewer->OnLigandsLoaded.RemoveAll(this);
    }

    // Bind to the OnLigandsLoaded event so we refresh when new molecules load
    CachedViewer->OnLigandsLoaded.AddDynamic(this, &UMoleculeListWidget::OnLigandsLoaded);

    // Populate immediately if there are already molecules
    Viewer->PopulateMoleculeListView(MoleculeListView);
}

void UMoleculeListWidget::AutoPopulate()
{
    // Find the first PDBViewer actor in the world
    APDBViewer* Viewer = Cast<APDBViewer>(
        UGameplayStatics::GetActorOfClass(GetWorld(), APDBViewer::StaticClass())
    );

    if (Viewer)
    {
        PopulateList(Viewer);
    }
}

void UMoleculeListWidget::OnLigandsLoaded()
{
    // Refresh the list when new ligands are loaded
    if (CachedViewer && MoleculeListView)
    {
        CachedViewer->PopulateMoleculeListView(MoleculeListView);
    }
}
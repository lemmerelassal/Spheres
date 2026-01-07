// PDBStructureWidget.cpp
#include "PDBStructureWidget.h"
#include "PDBViewer.h"
#include "Components/TreeView.h"
#include "Kismet/GameplayStatics.h"

void UPDBStructureWidget::NativeConstruct()
{
    Super::NativeConstruct();

    // Find the PDBViewer actor in the level
    TArray<AActor*> FoundActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundActors);
    
    if (FoundActors.Num() > 0)
    {
        PDBViewerRef = Cast<APDBViewer>(FoundActors[0]);
        
        if (PDBViewerRef)
        {
            // Bind to the OnResiduesLoaded event
            PDBViewerRef->OnResiduesLoaded.AddDynamic(this, &UPDBStructureWidget::OnStructureLoaded);
            
            // If TreeView is bound, set up delegates
            if (StructureTreeView)
            {
                // Bind the "get children" function - UE 5.6 syntax
                StructureTreeView->SetOnGetItemChildren(this, &UPDBStructureWidget::OnGetItemChildren);
                
                // Initial population (if structure already loaded)
                OnStructureLoaded();
            }
        }
    }
}

void UPDBStructureWidget::NativeDestruct()
{
    // Unbind from events
    if (PDBViewerRef)
    {
        PDBViewerRef->OnResiduesLoaded.RemoveDynamic(this, &UPDBStructureWidget::OnStructureLoaded);
    }
    
    Super::NativeDestruct();
}

void UPDBStructureWidget::OnStructureLoaded()
{
    if (!StructureTreeView || !PDBViewerRef)
        return;
    
    // Populate the tree view
    PDBViewerRef->PopulateTreeView(StructureTreeView);
}

void UPDBStructureWidget::OnGetItemChildren(UObject* Item, TArray<UObject*>& OutChildren)
{
    if (!PDBViewerRef)
        return;
    
    // Cast to our node type
    UPDBTreeNode* Node = Cast<UPDBTreeNode>(Item);
    if (!Node)
        return;
    
    // Get children from PDBViewer
    OutChildren = PDBViewerRef->GetChildrenForNode(Node);
}
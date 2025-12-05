#pragma once

#include "CoreMinimal.h"
#include "Widgets/SCompoundWidget.h"
#include "Widgets/DeclarativeSyntaxSupport.h"  // For using SLATE arguments

class ABlueSphere;

class FResidueVisibilityWidget : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(FResidueVisibilityWidget) {}
        SLATE_ARGUMENT(TWeakObjectPtr<ABlueSphere>, BlueSphere) // The BlueSphere actor
    SLATE_END_ARGS()

    void Construct(const FArguments& InArgs);

private:
    // A pointer to the BlueSphere instance
    TWeakObjectPtr<ABlueSphere> BlueSphere;

    // List to store buttons for each residue
    TArray<TSharedPtr<SButton>> Buttons;

    // Called when a residue's button is clicked to toggle visibility
    FReply OnResidueButtonClicked(int32 ResidueIndex);
};

#pragma once

#include "CoreMinimal.h"
#include "Widgets/SCompoundWidget.h"
#include "Widgets/DeclarativeSyntaxSupport.h"  // For using SLATE arguments

class ABlueSphere;

/**
 * Widget that displays a row of buttons for each residue.
 * Each button changes color depending on whether its residue is visible (blue) or hidden (red).
 */
class FResidueVisibilityWidget : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(FResidueVisibilityWidget) {}
        SLATE_ARGUMENT(TWeakObjectPtr<ABlueSphere>, BlueSphere) // The BlueSphere actor
    SLATE_END_ARGS()

    /** Builds the widget */
    void Construct(const FArguments& InArgs);

private:
    /** Reference to the BlueSphere actor */
    TWeakObjectPtr<ABlueSphere> BlueSphere;

    /** Buttons representing each residue */
    TArray<TSharedPtr<SButton>> ResidueButtons;

    /** Called when a residue's button is clicked to toggle visibility */
    FReply OnResidueButtonClicked(int32 ResidueIndex);

    /** Returns the button color depending on residue visibility (blue if shown, red if hidden) */
    FSlateColor GetResidueButtonColor(int32 ResidueIndex) const;
};

#pragma once

#include "CoreMinimal.h"
#include "Widgets/SCompoundWidget.h"
#include "Widgets/DeclarativeSyntaxSupport.h"  // For using SLATE arguments

class ALigands;

/**
 * Widget that displays a row of buttons for each Ligand.
 * Each button changes color depending on whether its Ligand is visible (blue) or hidden (red).
 */
class FLigandVisibilityWidget : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(FLigandVisibilityWidget) {}
        SLATE_ARGUMENT(TWeakObjectPtr<ALigands>, Ligands) // The Ligands actor
    SLATE_END_ARGS()

    /** Builds the widget */
    void Construct(const FArguments& InArgs);

private:
    /** Reference to the Ligands actor */
    TWeakObjectPtr<ALigands> Ligands;

    /** Buttons representing each Ligand */
    TArray<TSharedPtr<SButton>> LigandButtons;

    /** Called when a Ligand's button is clicked to toggle visibility */
    FReply OnLigandButtonClicked(int32 LigandIndex);

    /** Returns the button color depending on Ligand visibility (blue if shown, red if hidden) */
    FSlateColor GetLigandButtonColor(int32 LigandIndex) const;
};

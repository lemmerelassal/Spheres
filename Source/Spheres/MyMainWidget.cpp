// Fill out your copyright notice in the Description page of Project Settings.

#include "MyMainWidget.h"

void UMyMainWidget::GenerateButtons(int32 Count) {

    // Safety checks
    if (!ButtonContainer || !ButtonWidgetClass) {
        UE_LOG(LogTemp, Warning, TEXT("ButtonContainer or ButtonWidgetClass is null!"));
    }

    // Clear any existing buttons
    ButtonContainer->ClearChildren();

    for(int32 i = 0; i<Count; i++) {
        UUserWidget* NewButton = CreateWidget<UUserWidget>(this, ButtonWidgetClass);
        if(NewButton) {
            // Add button to conainer
            ButtonContainer->AddChild(NewButton);

            // Call Blueprint event for custom setup
            OnButtonCreated(NewButton, i);
        }
    }
}
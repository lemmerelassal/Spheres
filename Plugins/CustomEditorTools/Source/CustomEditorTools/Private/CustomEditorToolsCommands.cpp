// Copyright Epic Games, Inc. All Rights Reserved.

#include "CustomEditorToolsCommands.h"

#define LOCTEXT_NAMESPACE "FCustomEditorToolsModule"

void FCustomEditorToolsCommands::RegisterCommands()
{
	UI_COMMAND(PluginAction, "CustomEditorTools", "Execute CustomEditorTools action", EUserInterfaceActionType::Button, FInputChord());
	
	// Register the new commands
	UI_COMMAND(Button2Action, "Button 2 Action", "Execute action for button 2", EUserInterfaceActionType::Button, FInputChord());
	UI_COMMAND(Button3Action, "Button 3 Action", "Execute action for button 3", EUserInterfaceActionType::Button, FInputChord());
}

#undef LOCTEXT_NAMESPACE

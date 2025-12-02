// Copyright Epic Games, Inc. All Rights Reserved.

#include "CustomEditorTools.h"
#include "CustomEditorToolsStyle.h"
#include "CustomEditorToolsCommands.h"
#include "Misc/MessageDialog.h"
#include "ToolMenus.h"

static const FName CustomEditorToolsTabName("CustomEditorTools");

#define LOCTEXT_NAMESPACE "FCustomEditorToolsModule"

void FCustomEditorToolsModule::StartupModule()
{
	FCustomEditorToolsStyle::Initialize();
	FCustomEditorToolsStyle::ReloadTextures();

	FCustomEditorToolsCommands::Register();
	
	PluginCommands = MakeShareable(new FUICommandList);

	// Map the first button action
	PluginCommands->MapAction(
		FCustomEditorToolsCommands::Get().PluginAction,
		FExecuteAction::CreateRaw(this, &FCustomEditorToolsModule::PluginButtonClicked),
		FCanExecuteAction());

	// Map the second button action (NEW)
	PluginCommands->MapAction(
		FCustomEditorToolsCommands::Get().Button2Action,
		FExecuteAction::CreateRaw(this, &FCustomEditorToolsModule::Button2Clicked),
		FCanExecuteAction());

	// Map the third button action (NEW)
	PluginCommands->MapAction(
		FCustomEditorToolsCommands::Get().Button3Action,
		FExecuteAction::CreateRaw(this, &FCustomEditorToolsModule::Button3Clicked),
		FCanExecuteAction());

    // Store the handle when registering the callback (MODIFIED)
	StartupCallbackHandle = UToolMenus::RegisterStartupCallback(FSimpleMulticastDelegate::FDelegate::CreateRaw(this, &FCustomEditorToolsModule::RegisterMenus));
}
void FCustomEditorToolsModule::ShutdownModule()
{
	// This function may be called during shutdown to clean up your module.  For modules that support dynamic reloading,
	// we call this function before unloading the module.

    // REMOVE THIS LINE:
	// UToolMenus::UnregisterStartupCallback(StartupCallbackHandle);

	// KEEP THIS LINE:
	UToolMenus::UnregisterOwner(this);

	FCustomEditorToolsStyle::Shutdown();

	FCustomEditorToolsCommands::Unregister();
}


// Handler for the first button
void FCustomEditorToolsModule::PluginButtonClicked()
{
	FText DialogText = LOCTEXT("PluginButtonDialogText", "Button 1 clicked!");
	FMessageDialog::Open(EAppMsgType::Ok, DialogText);
}

// Handler for the second button (NEW)
void FCustomEditorToolsModule::Button2Clicked()
{
	FText DialogText = LOCTEXT("Button2DialogText", "Button 2 clicked!");
	FMessageDialog::Open(EAppMsgType::Ok, DialogText);
}

// Handler for the third button (NEW)
void FCustomEditorToolsModule::Button3Clicked()
{
	FText DialogText = LOCTEXT("Button3DialogText", "Button 3 clicked!");
	FMessageDialog::Open(EAppMsgType::Ok, DialogText);
}

void FCustomEditorToolsModule::RegisterMenus()
{
	// Owner will be used for cleanup in call to UToolMenus::UnregisterOwner
	FToolMenuOwnerScoped OwnerScoped(this);

	{
		UToolMenu* Menu = UToolMenus::Get()->ExtendMenu("LevelEditor.MainMenu.Window");
		{
			FToolMenuSection& Section = Menu->FindOrAddSection("WindowLayout");
			Section.AddMenuEntryWithCommandList(FCustomEditorToolsCommands::Get().PluginAction, PluginCommands);
		}
	}

	{
		UToolMenu* ToolbarMenu = UToolMenus::Get()->ExtendMenu("LevelEditor.LevelEditorToolBar.PlayToolBar");
		{
			FToolMenuSection& Section = ToolbarMenu->FindOrAddSection("PluginTools");
			
			// Add all three buttons to the toolbar section
			Section.AddEntry(FToolMenuEntry::InitToolBarButton(FCustomEditorToolsCommands::Get().PluginAction));
			Section.AddEntry(FToolMenuEntry::InitToolBarButton(FCustomEditorToolsCommands::Get().Button2Action));
			Section.AddEntry(FToolMenuEntry::InitToolBarButton(FCustomEditorToolsCommands::Get().Button3Action));
		}
	}
}

#undef LOCTEXT_NAMESPACE
	
IMPLEMENT_MODULE(FCustomEditorToolsModule, CustomEditorTools)


function [] = update_tab_context_menu(handleStruct)
    % UPDATE_TAB_CONTEXT_MENU - updates context menus associated with the
    %    inputted handle struct (for a fancy tab) as appropriate depending
    %    on the situation (e.g. what tabbed figures are available, how many
    %    siblings a tab currently has on each side, etc.) This is used as a
    %    callback for the main context menu so that options available are
    %    in sync with the active situation
    %
    % Inputs:
    %   handleStruct: a 1x1 struct containing various handles
    %      as values for its fields or as values for the fields in
    %      the substruct in the field shift. The struct can be generated by
    %      add_a_new_fancy_tab and contains the following information:
    %     hTab: the tab itself that was created
    %     hFig: the tabgroup's figure
    %     hFigContextMenu: the context menu on the
    %       figure
    %     hCloseTab: context menu option that allows one to close
    %       the tab
    %     hShiftContextMenu: context menu option that contains
    %       any suboptions to shift the tab left or right
    %     shift.hFarLeft: context menu option that allows one to
    %       shift the tab to the leftmost position in the tabgroup
    %     shift.hLeft: context menu option that allows one to
    %       shift the tab to the left in the tabgroup by one
    %       position
    %     shift.hRight: context menu option that allows one to
    %       shift the tab to the right in the tabgroup by one
    %       position
    %     shift.hFarLeft:  context menu option that allows one to
    %       shift the tab to the rightmost position in the tabgroup
    %     hRelocateContextMenu: context menu  that will contain
    %       any suboptions to relocate the tab
    %
    %  Side-effects:
    %    The handles provided, including that of the context menu, in the
    %      input struct are used to determine which options should be
    %      enabled/disabled/populated in the context menus based
    %      dynamically on the situation.
    %     Specifically  relocation options are generated dynamically
    %      based on what other tabbed figures are open at the time, and
    %      shifting options are enabled/disabled dynamically
    %      based on how many sibling tabs a tab has on each side
    %
    % Authors:
    %   Saair Quaderi
    
    import Fancy.UI.FancyTabs.delete_children;
    import Fancy.UI.FancyTabs.get_tab_shift_range;
    import Fancy.UI.FancyTabs.get_tab_relocation_options;

    if isfield(handleStruct, 'hRelocateContextMenu')
        [relocationCallbacks, relocationLabels] = get_tab_relocation_options(handleStruct.hTab);
        numRelocationOptions = numel(relocationLabels);
        delete_children(handleStruct.hRelocateContextMenu);
        for relocationOptionNum=1:numRelocationOptions
             uimenu(...
                'Parent', handleStruct.hRelocateContextMenu,...
                'Label', relocationLabels{relocationOptionNum},...
                'Callback', relocationCallbacks{relocationOptionNum});
        end
    end
    if isfield(handleStruct, 'hRelocateContextMenu')
        if numRelocationOptions > 0
            set(handleStruct.hRelocateContextMenu, 'Enable', 'on');
        else
            set(handleStruct.hRelocateContextMenu, 'Enable', 'off');
        end
    end
    [maxLeftShiftAmount, maxRightShiftAmount] = get_tab_shift_range(handleStruct.hTab);
    hasShiftOption = false;
    if isfield(handleStruct, 'shift')
        hMenuShift = handleStruct.shift;
        if maxLeftShiftAmount < 0
            if isfield(hMenuShift, 'hFarLeft')
                set(hMenuShift.hFarLeft, 'Enable', 'on');
            end
            if isfield(hMenuShift, 'hLeft')
                set(hMenuShift.hLeft, 'Enable', 'on');
            end
            hasShiftOption = true;
        else
            if isfield(hMenuShift, 'hFarLeft')
                set(hMenuShift.hFarLeft, 'Enable', 'off');
            end
            if isfield(hMenuShift, 'hLeft')
                set(hMenuShift.hLeft, 'Enable', 'off');
            end
        end
        if maxRightShiftAmount > 0
            if isfield(hMenuShift, 'hRight')
                set(hMenuShift.hRight, 'Enable', 'on');
            end
            if isfield(hMenuShift, 'hFarRight')
                set(hMenuShift.hFarRight, 'Enable', 'on');
            end
            hasShiftOption = true;
        else
            if isfield(hMenuShift, 'hRight')
                set(hMenuShift.hRight, 'Enable', 'off');
            end
            if isfield(hMenuShift, 'hFarRight')
                set(hMenuShift.hFarRight, 'Enable', 'off');
            end
        end
    end
    if isfield(handleStruct, 'hShiftContextMenu')
        if hasShiftOption
            set(handleStruct.hShiftContextMenu, 'Enable', 'on');
        else
            set(handleStruct.hShiftContextMenu, 'Enable', 'off');
        end
    end
end
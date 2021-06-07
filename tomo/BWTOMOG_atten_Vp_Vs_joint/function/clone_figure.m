%% Clone a figure
%  clone_figure()     - Clone the current figure to a new figure
%  clone_figure(N)    - Clone the current figure to figure N
%  clone_figure(0, M) - Clone figure M to a new figure
%  clone_figure(N, M) - Clone figure M to a figure N
function clone_figure(old_index, new_index)
    % Get handle of original figure
    if (exist('old_index', 'var'))
        if (not(ishandle(old_index)))
            error('myApp:argChk', ['Figure ',num2str(old_index),' undefined']);
        end
        figure(old_index);
    end
    old_handle = gcf;
    % Get handle of clone
    if (exist('new_index', 'var') && new_index)
        if (new_index < 0)
            error('myApp:argChk', ['Index ',num2str(new_index),' invalid']);
        end
        figure(new_index);
    else
        figure;
    end
    new_handle = gcf;
    % Copy figure's contents and properties
    copyobj(get(old_handle, 'Children'), new_handle);
    copyobj(get(old_handle, 'CurrentAxes'), new_handle);
    set(new_handle, 'Alphamap', get(old_handle, 'Alphamap'));
    set(new_handle, 'CloseRequestFcn', get(old_handle, 'CloseRequestFcn'));
    set(new_handle, 'Color', get(old_handle, 'Color'));
    set(new_handle, 'Colormap', get(old_handle, 'Colormap'));
    set(new_handle, 'CurrentPoint', get(old_handle, 'CurrentPoint'));
    set(new_handle, 'DockControls', get(old_handle, 'DockControls'));
    set(new_handle, 'FileName', get(old_handle, 'FileName'));
    set(new_handle, 'IntegerHandle', get(old_handle, 'IntegerHandle'));
    set(new_handle, 'InvertHardcopy', get(old_handle, 'InvertHardcopy'));
    set(new_handle, 'KeyPressFcn', get(old_handle, 'KeyPressFcn'));
    set(new_handle, 'KeyReleaseFcn', get(old_handle, 'KeyReleaseFcn'));
    set(new_handle, 'MenuBar', get(old_handle, 'MenuBar'));
    set(new_handle, 'Name', get(old_handle, 'Name'));
    set(new_handle, 'NextPlot', get(old_handle, 'NextPlot'));
    set(new_handle, 'NumberTitle', get(old_handle, 'NumberTitle'));
    set(new_handle, 'PaperUnits', get(old_handle, 'PaperUnits'));
    set(new_handle, 'PaperOrientation', get(old_handle, 'PaperOrientation'));
    set(new_handle, 'PaperPosition', get(old_handle, 'PaperPosition'));
    set(new_handle, 'PaperPositionMode', get(old_handle, 'PaperPositionMode'));
    set(new_handle, 'PaperSize', get(old_handle, 'PaperSize'));
    set(new_handle, 'PaperType', get(old_handle, 'PaperType'));
    set(new_handle, 'Pointer', get(old_handle, 'Pointer'));
    set(new_handle, 'PointerShapeCData', get(old_handle, 'PointerShapeCData'));
    set(new_handle, 'PointerShapeHotSpot', get(old_handle, 'PointerShapeHotSpot'));
    set(new_handle, 'Position', get(old_handle, 'Position'));
    set(new_handle, 'Renderer', get(old_handle, 'Renderer'));
    set(new_handle, 'RendererMode', get(old_handle, 'RendererMode'));
    set(new_handle, 'Resize', get(old_handle, 'Resize'));
    set(new_handle, 'ResizeFcn', get(old_handle, 'ResizeFcn'));
    set(new_handle, 'SelectionType', get(old_handle, 'SelectionType'));
    set(new_handle, 'ToolBar', get(old_handle, 'ToolBar'));
    set(new_handle, 'Units', get(old_handle, 'Units'));
    set(new_handle, 'WindowButtonDownFcn', get(old_handle, 'WindowButtonDownFcn'));
    set(new_handle, 'WindowButtonMotionFcn', get(old_handle, 'WindowButtonMotionFcn'));
    set(new_handle, 'WindowButtonUpFcn', get(old_handle, 'WindowButtonUpFcn'));
    set(new_handle, 'WindowKeyPressFcn', get(old_handle, 'WindowKeyPressFcn'));
    set(new_handle, 'WindowKeyReleaseFcn', get(old_handle, 'WindowKeyReleaseFcn'));
    set(new_handle, 'WindowScrollWheelFcn', get(old_handle, 'WindowScrollWheelFcn'));
    set(new_handle, 'WindowStyle', get(old_handle, 'WindowStyle'));
    set(new_handle, 'ButtonDownFcn', get(old_handle, 'ButtonDownFcn'));
    set(new_handle, 'Clipping', get(old_handle, 'Clipping'));
    set(new_handle, 'CreateFcn', get(old_handle, 'CreateFcn'));
    set(new_handle, 'DeleteFcn', get(old_handle, 'DeleteFcn'));
    set(new_handle, 'BusyAction', get(old_handle, 'BusyAction'));
    set(new_handle, 'HandleVisibility', get(old_handle, 'HandleVisibility'));
    set(new_handle, 'HitTest', get(old_handle, 'HitTest'));
    set(new_handle, 'Interruptible', get(old_handle, 'Interruptible'));
    set(new_handle, 'Parent', get(old_handle, 'Parent'));
    set(new_handle, 'Selected', get(old_handle, 'Selected'));
    set(new_handle, 'SelectionHighlight', get(old_handle, 'SelectionHighlight'));
    set(new_handle, 'Tag', get(old_handle, 'Tag'));
    set(new_handle, 'UIContextMenu', get(old_handle, 'UIContextMenu'));
    set(new_handle, 'UserData', get(old_handle, 'UserData'));
    set(new_handle, 'Visible', get(old_handle, 'Visible'));
    % from get(gcf) but skip CurrentCharacter, CurrentObject, BeingDeleted, Type
end
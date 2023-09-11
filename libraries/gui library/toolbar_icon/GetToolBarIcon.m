function cdata = GetToolBarIcon(icon)


switch icon
    
    case 'zoom in'
        
        cdata = imread('tool_zoom_in.png');

    case 'double tool arrow'
        
        cdata = imread('double_arrow.png');

    case 'pan'
        
        cdata = imread('tool_hand.png');
 
    case 'select node'
        
        cdata = imread('Edit_Node.png');
        
        
    case 'select element'
        
        cdata = imread('select_element_face.png');

    case 'bad element threshold'
        
        cdata = imread('bad_element.png');
        
    case 'zoom to bad element'
        
        cdata = imread('zoom_in_to_bad_element.png');
        
    case 'green arrow'
        
        cdata = imread('green_arrow.png');
        
    case 'create node string'

        cdata = imread('CreateNodeString.png');
        
    case 'google earth'
        
        [cdata,cmap] = imread('webicon.gif');
        
        cdata = ind2rgb(cdata,cmap);
        
        [i,j,~] = size(cdata);
        
        inan = find(cdata(:,:,1) == 1 & cdata(:,:,2) == 1 & cdata(:,:,3) == 1);
        
        BGFigColor = get(0,'DefaultUicontrolBackgroundColor');
        
        cdata(inan) = BGFigColor(1);
        cdata(inan+i*j) = BGFigColor(2);
        cdata(inan+i*j*2) = BGFigColor(3);
        return
        
    case 'pin'
        
        [cdata,cmap] = imread('pin_icon.gif');
        
        cdata = ind2rgb(cdata,cmap);
        
        [i,j,~] = size(cdata);
        
        inan = find(cdata(:,:,1) == 1 & cdata(:,:,2) == 1 & cdata(:,:,3) == 1);
        
        BGFigColor = get(0,'DefaultUicontrolBackgroundColor');
        
        cdata(inan) = BGFigColor(1);
        cdata(inan+i*j) = BGFigColor(2);
        cdata(inan+i*j*2) = BGFigColor(3);
        return
end

BGFigColor = get(0,'DefaultUicontrolBackgroundColor');

% Normalize BGFigColor to 0-255
newVal = ((BGFigColor - 0).*(255/(1-0)));

[i,j,~] = size(cdata);

inan = find(cdata(:,:,1) == 255 & cdata(:,:,2) == 255 & cdata(:,:,3) == 255);

cdata(inan) = newVal(1);
cdata(inan+i*j) = newVal(2);
cdata(inan+i*j*2) = newVal(3);
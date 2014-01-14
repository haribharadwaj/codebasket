classdef psychShape < handle
    % A shape class to be used to construct specific shapes for
    % visual cues such as arrows, fixation dots, response circles etc.
    % USAGE:
    %   SH = psychShape(PS,shapetype,shapesize,shapecentre)
    %       PS - The psychStarter object
    %       shapetype - Currently one of the following:
    %           'LARROW', 'RARROW', 'DOT', 'CIRCLE',
    %           'LRARROW', 'DOTR', 'DOTB', 'DOTP'
    %       shapesize - In pixel units (see properties of PS)
    %       shapecentre - In pixel units (see properties of PS)
    %
    %--------------------------
    % Hari Bharadwaj, September 6, 2010
    %--------------------------
    properties (SetAccess = private)
        shape % The actual type of shape it is
        size % Size of the bounding box (in pixel coordinates)
        centre % Location of the bounding box (in pixel coordinates): 2D
        renderArgs % Arguments needed to render the shape (cells)
        command % Command for rendering the shape usually called with Screen();
        colour
        thickness
        intensityFactor
    end
    
    events
        PTBbadVersion % Relies on PTB textures appropriate version needed
    end
    
    methods
        function SH = psychShape(PS,shapetype,shapesize,shapecentre)
            SH.shape = shapetype;
            SH.size = shapesize;
            SH.centre = shapecentre;
            SH.colour = PS.white;
            SH.intensityFactor = 0.5;
            
            switch(SH.shape)
                case 'LARROW'
                    SH.thickness = round(SH.size/6);
                    SH.renderArgs{1} =  +...
                        [[0; SH.size]...
                        [-SH.size; 0]...
                        [-SH.size; 0]...
                        [0; -SH.size]];
                    SH.command = 'DrawLines';
                    SH.renderArgs{2} = SH.thickness;
                    SH.renderArgs{3} = SH.colour*SH.intensityFactor;
                    SH.renderArgs{4} = SH.centre;
                case 'RARROW'
                    SH.thickness = round(SH.size/6);
                    SH.renderArgs{1} =  +...
                        [[0; SH.size]...
                        [SH.size; 0]...
                        [SH.size; 0]...
                        [0; -SH.size]];
                    SH.command = 'DrawLines';
                    SH.renderArgs{2} = SH.thickness;
                    SH.renderArgs{3} = SH.colour*SH.intensityFactor;
                    SH.renderArgs{4} = SH.centre;
                case 'DOTB'
                    SH.thickness = SH.size;
                    SH.colour = PS.blue*SH.intensityFactor;
                    SH.renderArgs{1} = SH.colour;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                case 'DOTR'
                    SH.thickness = SH.size;
                    SH.colour = PS.red*SH.intensityFactor;
                    SH.renderArgs{1} = SH.colour;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                case 'DOTP'
                    SH.thickness = SH.size;
                    SH.colour = 0.5*(PS.red+PS.blue)*SH.intensityFactor;
                    SH.renderArgs{1} = SH.colour;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                case 'DOT'
                    SH.thickness = SH.size;
                    SH.renderArgs{1} = SH.colour*SH.intensityFactor;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                case 'DOTX'
                    SH.thickness = SH.size;
                    SH.renderArgs{1} = SH.colour;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                    
                case 'CIRCLE'
                    SH.thickness = round(SH.size/8);
                    SH.renderArgs{1} = SH.colour*SH.intensityFactor;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                    
                case 'GCIRCLE'
                    SH.thickness = round(SH.size/8);
                    SH.renderArgs{1} = PS.green*SH.intensityFactor;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                case 'BCIRCLE'
                    SH.thickness = round(SH.size/8);
                    SH.renderArgs{1} = PS.blue*SH.intensityFactor;
                    SH.renderArgs{2} = SetRect(...
                        SH.centre(1) - SH.size,... % Left
                        SH.centre(2) - SH.size,... % Top
                        SH.centre(1) + SH.size,... % Right
                        SH.centre(2) + SH.size); % Bottom
                    SH.command = 'FrameOval';
                    SH.renderArgs{3} = SH.thickness;
                    
                case 'RCROSS'
                    SH.thickness = round(SH.size/5);
                    SH.renderArgs{1} =  +...
                        [[-SH.size; SH.size]...
                        [SH.size; -SH.size]...
                        [-SH.size; -SH.size]...
                        [SH.size; SH.size]];
                    SH.command = 'DrawLines';
                    SH.renderArgs{2} = SH.thickness;
                    SH.renderArgs{3} = PS.red*SH.intensityFactor;
                    SH.renderArgs{4} = SH.centre;
                    
                case 'LRARROW'
                    SH.thickness = round(SH.size/6);
                    SH.renderArgs{1} =  +...
                        [[0; SH.size]...
                        [-SH.size; 0]...
                        [-SH.size; 0]...
                        [0; -SH.size]...
                        [0; SH.size]...
                        [SH.size; 0]...
                        [SH.size; 0]...
                        [0; -SH.size]];
                    SH.command = 'DrawLines';
                    SH.renderArgs{2} = SH.thickness;
                    SH.renderArgs{3} = SH.colour*SH.intensityFactor;
                    SH.renderArgs{4} = SH.centre;
                otherwise
                    fprintf(1,'\nERROR! Unkown Shape Type! Cannot create shape.. Sorry!\n');
            end
            
        end
        function renderShape(SH,PS)
            switch(numel(SH.renderArgs))
                case 1
                    Screen(SH.command,PS.window,SH.renderArgs{1});
                case 2
                    Screen(SH.command,PS.window,SH.renderArgs{1},...
                        SH.renderArgs{2});
                case 3
                    Screen(SH.command,PS.window,SH.renderArgs{1},...
                        SH.renderArgs{2},SH.renderArgs{3});
                case 4
                    Screen(SH.command,PS.window,SH.renderArgs{1},...
                        SH.renderArgs{2},SH.renderArgs{3}, SH.renderArgs{4});
                otherwise
                    sca;
                    fprintf(1,'\nCould not render shape! Sorry!\n');
            end
        end
    end
    
end

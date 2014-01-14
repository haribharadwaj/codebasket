function renderVisFrame(PS,frametype)
% A function to generate different frames such as the
% cue frame, fixation frame etc.
% It mostly aggregates different psychShapes. The frame is rendered
% when it is being created.
% USAGE:
%   VFR = psychVisFrame(PS, frametype)
%       PS - The psychStarter object
%       frametype - Currently one of the following:
%           'CUEL', 'CUER', 'FIX', 'RESP', 'CUEC', 'CORRECT', 'WRONG',
%           'NORESP'
%
%--------------------------
% Hari Bharadwaj, September 6, 2010
%--------------------------


visangle = 1;
dotangle = 0.15;

switch(frametype)
    case 'CUEL'
        LA = psychShape(PS,'LARROW',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        LA.renderShape(PS);
        DOT.renderShape(PS);
    case 'CUER'
        RA = psychShape(PS,'RARROW',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        RA.renderShape(PS);
        DOT.renderShape(PS);
    case 'CUEC'
        LRA = psychShape(PS,'LRARROW',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        LRA.renderShape(PS);
        DOT.renderShape(PS);
    case 'FIX'
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        DOT.renderShape(PS);
    case 'FLICKER'
        DOT = psychShape(PS,'DOTX',dotangle*PS.oneDeg,PS.centre);
        DOT.renderShape(PS);
    case 'CORRECT'
        DOT = psychShape(PS,'DOTB',dotangle*PS.oneDeg,PS.centre);
        DOT.renderShape(PS);
    case 'WRONG'
        DOT = psychShape(PS,'DOTR',dotangle*PS.oneDeg,PS.centre);
        DOT.renderShape(PS);
    case 'NORESP'
        DOT = psychShape(PS,'DOTP',dotangle*PS.oneDeg,PS.centre);
        DOT.renderShape(PS);
    case 'RESP'
        CIRC = psychShape(PS,'CIRCLE',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        CIRC.renderShape(PS);
        DOT.renderShape(PS);
    case 'GO'
        GCIRC = psychShape(PS,'GCIRCLE',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        GCIRC.renderShape(PS);
        DOT.renderShape(PS);
    case 'GOBLUE'
        GCIRC = psychShape(PS,'BCIRCLE',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        GCIRC.renderShape(PS);
        DOT.renderShape(PS);
    case 'NOGO'
        RCROSS = psychShape(PS,'RCROSS',visangle*PS.oneDeg,PS.centre);
        DOT = psychShape(PS,'DOT',dotangle*PS.oneDeg,PS.centre);
        RCROSS.renderShape(PS);
        DOT.renderShape(PS);
        
end

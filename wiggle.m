function [handle,scale] = wiggle(seismic,time,offset,sc,ornt)
%WIGGLE				[ A.K.Booer 24-Mar-1993 ]
%
%   wiggle(seismic,time,offset) - plots a shaded wiggle plot
%
%   seismic - matrix of waveform columns
%   time - time axis vector
%   offset - space axis vector
% 
%   wiggle(seismic,t,z,mag)   - scales waveforms by given magnification
%   wiggle(seismic,t,z,m,'a') - changes time axis orientation to 'across'
%
%   WIGGLE returns handle to created graphics objects, so that
%   set(wiggle(x,t,z), 'face',c1, 'edge',c2') will specify colours,
%   and scale factor used

% Modifications:
% 24-Mar-1993   akb   Based on "wp" of 2-Oct-1992
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%seismic=seismic(:,1:10:end)
[n,m] = size(seismic);
n1 = n:-1:1;
if nargin < 2, time = 1:n; end
if nargin < 3, offset = 1:m; end
if nargin < 4, sc = 1; end
if nargin < 5, ornt = 'd' ; end		% default orientation (d)own/(a)cross
colour = 'k'; 				% default colour
offset = offset(:)';
% 
scale = (2*mean(diff(offset))) * (sc / max(max(seismic) - min(seismic)));
if sc < 0, scale= -sc; end;
seismic = seismic * scale;
%

if ornt == 'd'
  h = fill(offset(ones(2*n,1),:)+[seismic;min(seismic(n1,:),0)],time([1:n n1]),colour);
  set(h,'EdgeColor',colour), 
   set(gca,'Ydir','r')
  dz = offset(2)-offset(1);
  axis([min(offset)-dz max(offset)+dz min(time) max(time)]);
else
  h = fill(time([1:n n1]),offset(ones(2*n,1),:)-[seismic;min(seismic(n1,:),0)],colour);...
      set(h,'EdgeColor',colour), set(gca,'Ydir','r')
end


%
if nargout > 0, handle = h; end
%


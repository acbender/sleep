function [Rho, Pval, slope, offset] = regresscirclin (phases, position)

% circular linear regression of SO phases of spindles vs EEG channel position

% modified from Kempter, Schmidt 2012

% requires all of the following inputs:
% SO phases of spindle events (in radians)
% position vector of integers, corresponding to EEG channel positions along anterior-posterior axis (1=most anterior)

% outputs: 
% R, P, slope (in radians per x-unit), offset (y-intercept)  


%% check
    if length(phases)-length(position)~=0
        error('error: length of inputs must be equal')
    end

%% 

% Get rid of all the nans in this data and relabel
    circ=phases;
    lin=position;
    circ = circ(~isnan(phases));
    lin = lin(~isnan(phases));

% Make sure there are still valid data 
    if length(lin)==0 || length(circ)==0
        Rho=NaN; Pval=NaN; slope=NaN; offset=NaN; b=NaN;
        return 
    end

% finding the optimal slope, note that we have to restrict the range of slopes 
    slope_bounds=[-3*pi 0]; % slope range: must be restricted for optimization, [-3*pi 0] works best
    minslope=slope_bounds(1) / (max(lin)-min(lin));
    maxslope=slope_bounds(2) / (max(lin)-min(lin));
    p=[minslope:abs(maxslope-minslope)/1000:maxslope]';

    R=NaN([length(p) 1]);
    for i=1:length(p)
    R(i)=sqrt( (sum(cos(circ-(p(i)*lin)))/length(circ))^2 + (sum(sin(circ-(p(i)*lin)))/length(circ))^2 );
    end
    
    [m, m_i]=max(R);
    slope=p(m_i);

% calculate offset
    offset = atan2(sum(sin(circ-(slope*lin))), sum(cos(circ-(slope*lin))));  

% circular-linear correlation:
    [Rho, Pval] = circ_corrcl(circ, lin); 

    % Assign the correct sign to rho
    if slope < 0
        Rho = -abs(Rho);
    else
        Rho = abs(Rho);
    end


% END

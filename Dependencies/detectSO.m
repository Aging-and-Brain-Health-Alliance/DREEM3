% %----------------------------------------------------------------------%
% function: detectSO
%
%           Function to detect slow oscillations, similar to a range of 
%           published studies, e.g. Massimini et al. 2004. First the raw
%           signal is filtered in a pre-set range (0.1 - 4 Hz, butterworth
%           filters, low-pass of order 10, high-pass of order 4, using
%           filtfilt for forward and backward pass so as to not distort the
%           phase). Then, all positive-to-negative zero crossings are 
%           detected (save for a buffering 2 sec at the beginning and end
%           of recording) and the following half-wave until the next
%           negative-to-positive zero crossing is a candidate event.
%           Finally, events are screened to satisfy a range of basic 
%           criteria - depth of slow oscillation, duration of trough, etc.
%
%           For example, to replicate Massimini et al. 2004:
%           detectSO(signal,freq,80,140,0.3,1,x) 
%           where x is large enough to always be satisfied (they did not
%           use the criterion of total length)
%
% dependancy: - filtering parameters and functions (see the beginning of
%               the code)
%
% input:   - raw EEG signal [uV];
%          - frequency of EEG acquisition [Hz];
%          - threshold depth of the voltage minimum [uV];
%          - threshold height of the voltage peak-to-peak, i.e. depth of of
%          the negative half-wave plus height of positive half-wave [uV];
%          - minimum duration of the negative half-wave [sec];
%          - maximum duration of the negative half-wave [sec];
%          - maximum duration of the whole wave [sec].
%
% output:   - matrix holding the timepoints of the beginning, middle and
%           end of the slow oscillation (i.e. appropriate zero crossings)
%           [in samples], each row is a SO, with respective time points in
%           columns;
%           - array holding the timepoints of the minimum trough of each 
%           detected slow oscillation (minimum of raw signal) [samples].
%
% WARNING: detection is on the filtered signal, therefore returned zero 
%          crossings might not correspond exactly to the raw signal zero
%          crossings (the filtering is to smooth out low-amplitude sharp
%          fluctuations that can move the zero crossings and thus distort
%          the width of the trough, which sometimes results in false
%          negatives)
%
%           written by Mohsen Naji, modified by Negin Sattari
% %----------------------------------------------------------------------%

function [can_events,events_minima]=detectSO(signal_raw,freq,threshold_neg,peak_to_peak,duration_min,duration_max,duration_whole)

% here pre-processing of EEG, filtering between 0.1-4 Hz:
% design the butterworth filters to use:
[b1,a1] = butter(4,0.1./(freq/2),'high');
% [b2,a2] = butter(10,4./(freq/2),'low');
[b2,a2] = butter(4,4./(freq/2),'low'); % Maryam
% and filter the raw signal:
x = filtfilt(b1,a1,signal_raw);
signal = filtfilt(b2,a2,x);


% initialize output as empty:
can_events = [];
events_minima = [];

zero_cross=[];
% go through the filtered signal to find all zero crossings:
for i = 2*freq+1:length(signal)-2*freq
    % easy way to check for zero crossins:
    if  signal(i)*signal(i+1)<0
        zero_cross=[zero_cross i+1];
    end
end

if length(zero_cross)>1
    
    % make sure the first zero crossing is positive-to-negative:
    if signal(zero_cross(1))<signal(zero_cross(1)+1)
       zero_cross(1)=[];
    end
    
    % cut down to the number of 'whole' candidate segments (negative
    % half-wave and a following positive half-wave):
    noc=floor(length(zero_cross)/3); 
    zero_cross=zero_cross(1:3*noc);
    % extract the putative beginning, middle and end:
    c1=zero_cross(1:2:end)';   % positive-to-negative zero crossing beginning the slow oscillation
    c2=zero_cross(2:2:end)';   % negative-to-positive zero crossing as the 'middle', end of negative half-wave
    c3=zero_cross(3:2:end)';   % positive-to-negative zero crossing ending the slow oscillation
    l=min([length(c1) length(c2) length(c3)]);
    % an array with the candidate segments:
    can_events=[c1(1:l) c2(1:l) c3(1:l)];
    
    % auxilliary variable with widths of the negative and positive
    % half-waves:
    prop=[can_events(:,2)-can_events(:,1) can_events(:,3)-can_events(:,2)];
    % select only the candidate events satisfying desired width criteria:
    can_events = can_events( prop(:,1)>duration_min*freq & prop(:,1)<duration_max*freq & prop(:,2)<duration_whole*freq ,:);
    clear prop
    
    if ~isempty(can_events)
        
        [rr,~] = size(can_events);
        prop = zeros(rr,2);
        % go through all candidate events and extract the minimum voltage
        % of the negative half-wave and maximum voltage of the positive
        % half-wave:
        for i=1:rr
            prop(i,1) = abs(min(signal(can_events(i,1):can_events(i,2))));
            prop(i,2) = abs(max(signal(can_events(i,2):can_events(i,3))));
        end
        
    end
    
    if ~isempty(can_events)
        % select only the candidate events satisfying desired voltage value
        % criteria:
        can_events = can_events( prop(:,1)>threshold_neg & (prop(:,1)+prop(:,2))>peak_to_peak ,:);
    end
    
    if ~isempty(can_events)
        
        [rr,~] = size(can_events);
        events_minima = zeros(rr,1);
        % go through all detected events and extract the position and value
        % of minimal RAW voltage (as opposed to filtered!):
        for i=1:rr
            ev_mini = min(signal_raw(can_events(i,1):can_events(i,2)));
            pos = find(signal_raw(can_events(i,1):can_events(i,2))==ev_mini);
            events_minima(i) = can_events(i,1) + pos(1);
        end
        
    end
    
end

end

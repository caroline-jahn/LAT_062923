function [Coupled_channels_with_sig_responce] = add_coupled_channels(ActiveChannels,Channels_with_sig_response)
%Takes the channels with sign activation to xx and add the coupled channel
%if it is also in the active channels list

if isempty(Channels_with_sig_response)==1
Coupled_channels_with_sig_responce=NaN;
end
    

count=1;
for i=1:length(ActiveChannels)
    if ismember(ActiveChannels(i),Channels_with_sig_response)==1
        Coupled_channels_with_sig_responce(count)=ActiveChannels(i); %if sig response add
        count=count+1;
    elseif mod(ActiveChannels(i),2)==0 && ismember(ActiveChannels(i)-1,Channels_with_sig_response)==1 %even channel (n) => look if coupled (n-1) is sig
        Coupled_channels_with_sig_responce(count)=ActiveChannels(i); %if sig response add
        count=count+1;
    elseif mod(ActiveChannels(i),2)==1 && ismember(ActiveChannels(i)+1,Channels_with_sig_response)==1 %odd channel (n) => look if coupled (n+1) is sig
        Coupled_channels_with_sig_responce(count)=ActiveChannels(i); %if sig response add
        count=count+1;
    end
end


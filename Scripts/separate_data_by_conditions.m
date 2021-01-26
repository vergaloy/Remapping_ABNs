function [mice_sleep,N]=separate_data_by_conditions(data)
N=[];
for i=1:size(data,1)
  temp=full(data{i, 1});
  N=[N;ones(size(temp,1),1)*i];
  mice_sleep(i,:)=Separate_data(temp,50);
end

end

function C=Separate_data(obj,remove_first_points)  % some time artifacts appear in the begining of a video sessions
%Because of this we remove the first 10 seconds
% C=Separate_by_behaviour_no_bin(neuron.S,hypno,50);

i=0;

temp=full(obj);
i=i+1;
C{i}=temp(:,1+remove_first_points:3000); %HC
i=i+1;
C{i}=temp(:,3001+remove_first_points:6000); %preS
i=i+1;
C{i}=temp(:,6001+remove_first_points:7500);  %postS
i=i+1;
C{i}=temp(:,7501+remove_first_points:52500); %Consolidaiton
i=i+1;
C{i}=temp(:,52501+remove_first_points:55500); %Test
end

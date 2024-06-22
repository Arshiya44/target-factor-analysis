
%%% TFA for finding k %%%%
%%%% HA -> H+ + A- %%%%
nuke
C_HA = 2;
k = 5*10^-5;
PH = 0:.2:14;
n = 1;
A = ones([1,length(PH)]);
HA = ones([1,length(PH)]);
for i = PH
    A(n)= (k*C_HA)/((10^-i)*C_HA+k);
    HA(n) = C_HA - A(n);
    n = n+1;
end
concentration_profiles = [HA;A];
figure('Name','concentration profuiles')
plot(PH,[HA;A])
xlabel('PH') ; ylabel('concentration')
name =["HA","A"];
legend(name)


base_line = 400:1200;

s_peaks = [800,900];
s_widths = [70,60];
s_scalers = [.7,.6];
pure_spectrum = ones([length(s_peaks),length(base_line)]);
for i = 1:length(s_peaks)
    pure_spectrum(i,:) = s_scalers(i)*gaussmf(base_line,[s_widths(i),s_peaks(i)]);
end

spectrum_data = concentration_profiles' * pure_spectrum;
[u, s, v]= svd(spectrum_data);


for i = 1:2
    vs(i,:) = v(:,i)*s(i,i);
end
figure('Name','u space')
scatter(vs(1,:),vs(2,:))
hold on


%%% simulated targets %%%
m =1;
k_range = 0:.2:7;
for ii = k_range
    n =1;
    k =  5*10^-ii;
    for i = PH
        sim_A(m,n)= (k*C_HA)/((10^-i)*C_HA+k);
        n = n+1;
    end
    m = m+1;
end
figure('Name','simulated cooncentration profiles')
plot(PH,sim_A)
hold on


%%

figure(2)
true_HA = HA * u(:,1:2) ;
true_A = A * u(:,1:2);
true = [true_HA;true_A];
scatter(true(:,1),true(:,2),300,'red','hexagram','LineWidth',1)
text_true = text(true(:,1),true(:,2),name,"HorizontalAlignment","center","FontSize",9);
sim_projection = sim_A * u(:,1:2);
scatter(sim_projection(:,1),sim_projection(:,2),'*')


figure('Name','reconstructed concentration profiles')
reconstructed = sim_projection * u(:,1:2)';
plot(PH,reconstructed)
n = 1;
for i = 1:length(k_range)
    residual = sim_A(i,:) - reconstructed(i,:);
    RSD(n) = sqrt(sum(svd(residual))^2/(2*length(k_range)));
    n = n+1;
end
figure('Name','rsd')
plot(k_range,RSD,'LineWidth',1)
xlabel('10^-')




















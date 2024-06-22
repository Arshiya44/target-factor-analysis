%%%% target factor analysis %%%%
nuke
%%%% concentration profiles %%%%%%%

time = 1:360;
c_peaks = [60,120];
c_widths = [22,32];
c_scalers = [.6,.4];
for i = 1:length(c_peaks)
    concentration_profiles(i,:) = c_scalers(i)*gaussmf(time,[c_widths(i),c_peaks(i)]);
end

%%%% spectrum profiles %%%%%%%%
base_line = 400:1200;

s_peaks = [800,900];
s_widths = [70,60];
s_scalers = [.7,.6];

for i = 1:length(s_peaks)
    pure_spectrum(i,:) = s_scalers(i)*gaussmf(base_line,[s_widths(i),s_peaks(i)]);
end


%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%

hplcdad = concentration_profiles' * pure_spectrum;
figure('Name','pure components')
subplot(2,1,1)
plot(base_line,pure_spectrum)
title('spectrum')
subplot(2,1,2)
plot(time,concentration_profiles)
title('concentration')
[loadings, scores] = pca(hplcdad,"Algorithm","svd","Centered",false,"Economy",false,"NumComponents",2);
figure('Name','v sapce')
scatter(scores(:,1),scores(:,2),'magenta','.')
hold on

%%% 8 pure spectrum %%%
base_line = 400:1200;

peaks = [s_peaks(1),s_peaks(2),700,730,900,1000,500,1100];
widths = [s_widths(1),s_widths(2),30,30,30,40,15,17];
scalers = [s_scalers(1),s_scalers(2),.5,.6,.7,.8,.2,.4];
spectrum_box = zeros([length(peaks),length(base_line)]);
for i = 1:length(peaks)
    spectrum_box(i,:) = scalers(i)*gaussmf(base_line,[widths(i),peaks(i)]);
end
figure('Name','spectrum box')
plot(base_line,spectrum_box(1:2,:),'LineStyle','-.','LineWidth',1)
hold on
plot( base_line,spectrum_box(3:end,:),'LineStyle','-','LineWidth',1)
name = strsplit(num2str(1:length(peaks)));
legend(name)


%%% projection of spectruns into v space %%%

projection_in_v_space = spectrum_box * loadings;
figure(2)
colors = rand(length(peaks), 3);
scatter(projection_in_v_space(:,1),projection_in_v_space(:,2),100,colors,"LineWidth",2);
text1 = text(projection_in_v_space(:,1),projection_in_v_space(:,2),name,"HorizontalAlignment","center","FontSize",8);

%%
reconstructed_spectrums = projection_in_v_space * loadings';
figure('Name','8 pure spectrum')
hold_size = round(length(peaks)/2);
for i = 1:length(peaks)
    subplot(hold_size,2,i)
    plot(base_line,spectrum_box(i,:),"Marker","*")
    hold on
    plot(base_line,reconstructed_spectrums(i,:),"Marker","v")
    title(num2str(i))
    hold_corelation = correlationvec(spectrum_box(i,:), ...
        reconstructed_spectrums(i,:));
    subtitle(['corelation = ',num2str(hold_corelation)])
end

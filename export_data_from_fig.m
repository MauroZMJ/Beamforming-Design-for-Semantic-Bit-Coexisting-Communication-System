clc;
clear all;
close all;
syms x
%h=open('./figures/Pareto_semantic_SNR_5.fig');
%h=open('./figures/Pareto_digital_SNR_5.fig');
%h=open('./figures/beamforming_comparison_differ_true_snr.fig');
h=open('./figures/k_optimization_case_2.fig');

lh=findall(gca,'type','line');
xc=get(lh,'xdata');
yc=get(lh,'ydata');


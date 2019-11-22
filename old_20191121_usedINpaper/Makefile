clean:clean_exe_read_chi2true_chi2fit clean_exe_plot_FC_Coverage_Interval clean_exe_read_FC_CL
	-rm -f *.d *.so *~

clean_exe_read_chi2true_chi2fit:
	-rm exe_read_chi2true_chi2fit	
exe_read_chi2true_chi2fit:clean_exe_read_chi2true_chi2fit
	g++ -o  exe_read_chi2true_chi2fit exe_read_chi2true_chi2fit.cc -O `root-config --libs` -I. -I/home/xji/software/root/include -L/home/xji/software/root/lib -lMathMore

clean_exe_plot_FC_Coverage_Interval:
	-rm exe_plot_FC_Coverage_Interval	
exe_plot_FC_Coverage_Interval:clean_exe_plot_FC_Coverage_Interval
	g++ -o  exe_plot_FC_Coverage_Interval exe_plot_FC_Coverage_Interval.cc -O `root-config --libs` -I. -I/home/xji/software/root/include -L/home/xji/software/root/lib -lMathMore

clean_exe_read_FC_CL:
	-rm exe_read_FC_CL	
exe_read_FC_CL:clean_exe_read_FC_CL
	g++ -o  exe_read_FC_CL exe_read_FC_CL.cc -O `root-config --libs` -I. -I/home/xji/software/root/include -L/home/xji/software/root/lib -lMathMore

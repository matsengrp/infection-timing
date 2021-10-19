import TI_predictor as predict
import pandas as pd

infant_datapath = '../_ignore/neher_webtool_results.tsv'
infant_data = pd.read_csv(infant_datapath, sep='\t')

apd_list = infant_data['average_APD'].tolist()
regions = list(zip(map(int, infant_data['HXB2nt_start']), map(int, infant_data['HXB2nt_end '])))

estimated_times = []
for i in range(len(apd_list)):
    time = predict.TI_from_diversity(apd_list[i], regions[i], 0.01)
    estimated_times.append(time[0])

infant_data['adult_model_time_since_infection_estimates'] = estimated_times

infant_data.to_csv(infant_datapath, sep='\t')

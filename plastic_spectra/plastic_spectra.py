import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks as fp

# This module creates database with 3 tables. One for plastic names, one for plastic residue, one for experemental spectras
def create_database():
    con = sqlite3.connect("plastic_spectrum.db")
    cur = con.cursor()
    try:
        cur.execute("""CREATE TABLE plastic(
                plastic_id INTEGER PRIMARY KEY,
                plastic_name VARCHAR(70)
                );""")
    except:
        print('Error plastic')
    try:
        cur.execute(""" CREATE TABLE residue(
                residue_id INTEGER PRIMARY KEY,
                mz DECIMAL(12,6),
                name VARCHAR(200),
                plastic_id INT,
                decoy_status INT,
                FOREIGN KEY (plastic_id)  REFERENCES plastic (plastic_id) ON DELETE CASCADE);""")
    except:
        print('Error residue')
    try:
        cur.execute("""CREATE TABLE spectra(
                peak_id INTEGER PRIMARY KEY,
                residue_id INT,
                name VARCHAR(200),
                mz DECIMAL(12,6),
                intensity DECIMAL(12,6),
                status VARCHAR(200));""")
        print('database plastic_spectrum.db successfully created')
    except:
        print('Error spectra')
    cur.close()
    con.close()




class plastic_spectra:


    def __init__(self):
        self.annotation = pd.DataFrame(columns = ["mz", "intensity", "residue_id", "status", "plastic_name", "description","match_mz"])

    def read_spectra(self, filename = 'any.txt', retention_time = 0.111111, threshold = 0.5):
        """This fuction reads spectras from ms1 format files,
        filename - any.txt with experemental results,
        retention_time - retention time of spectra in file"""
        with open(filename) as f: 
            data = f.readlines()
            data = [row.replace('\n', '') for row in data]
            scan_id = int(data[data.index(f'I\tRTime\t{retention_time}') - 1][data[data.index(f'I\tRTime\t{retention_time}') - 1].find('=') + 1:])
            data = data[data.index(f'I\tRTime\t{retention_time}') + 4 : data.index(f'S\t{scan_id + 1}\t{scan_id + 1}')]
            data = [i.split(' ') for i in data]
            data = [[float(j) for j in i] for i in data]
            data = pd.DataFrame(data, columns = ['mz', 'intensity'])
            filtered_data = pd.DataFrame(index = data.index)
            for i, vals in data.items():
                peaks, _ = fp(vals, threshold = threshold)
                filtered_data[i+'_peaks'] = filtered_data.index.isin(peaks + data.index.min())
                peak_filter = filtered_data.sum(axis = 1) >= 1
            self.annotation[["mz","intensity"]] = data.loc[peak_filter]


    def __residue_search(self, mz, tolerance):
        # Search in database for matching plastics based on m/z and tolerance gap
        con = sqlite3.connect("plastic_spectrum.db")
        cur = con.cursor()
        try:
            query = cur.execute("""SELECT re.residue_id as residue_id, pl.plastic_name as plastic_name,
                                re.name as re_name, re.mz as match_mz, re.decoy_status as decoy_status
                                FROM residue AS re JOIN
                                plastic AS pl on re.plastic_id = pl.plastic_id
                                WHERE re.mz BETWEEN ?-? AND ?+?""", (mz, tolerance, mz, tolerance,))
            cols = [column[0] for column in query.description]
            result = query.fetchall()
            return [{col : prop for col, prop in zip(cols, plastic)} for plastic in result]
        except:
            print('Database connection error')


    def __clean(self, m, x):
        # If don't get any results for given m/z in our database
        if x == []:
            return {'residue_id': None,
                    'plastic_name': None,
                    're_name': None,
                    'match_mz': None,
                    'decoy_status': None,
                    'status' : 'No matches'}
        # If all matches are decoys
        elif all(i['decoy_status']==1 for i in x):
            values = [abs(m-i['match_mz']) for i in x]
            result = x[values.index(min(values))]
            result['status'] = 'Random match'
            return result
        # If any match is decoy
        elif any(i['decoy_status']==1 for i in x):
            decoy_plastics = set([i['plastic_name'] for i in x if i['decoy_status'] == 1])
            not_decoy = [i for i in x if i['decoy_status'] == 0]
            without_decoy = [i for i in not_decoy if i['plastic_name'] not in decoy_plastics]
            if len(without_decoy) !=0:
                values = [abs(m-i['match_mz']) for i in without_decoy]
                result = without_decoy[values.index(min(values))]
                result['status'] = 'Match'
                return result   
            else:
                values = [abs(m-i['match_mz']) for i in not_decoy]
                result = not_decoy[values.index(min(values))]
                result['status'] = 'Possibility of random match'
                return result   
        # If there is NO decoys
        else:
            values = [abs(m-i['match_mz']) for i in x]
            result = x[values.index(min(values))]
            result['status'] = 'Match'
            return result
        

    def annotate_spectra(self, tolerance = 0.1):
        """ Annotates your spectra based on tolerance factor"""
        try:
            row_annotation = {i : self.__residue_search(i, tolerance) for i in self.annotation['mz']} # Get all matching plastics from database
            clean_annotation = {i : self.__clean(i, row_annotation[i]) for i in row_annotation.keys()} # Choose closest plastic, and adds status of match
            self.annotation["residue_id"] = self.annotation["mz"].apply(lambda x: clean_annotation[x]['residue_id'])
            self.annotation["status"] = self.annotation["mz"].apply(lambda x: clean_annotation[x]['status'])
            self.annotation["plastic_name"] = self.annotation["mz"].apply(lambda x: clean_annotation[x]['plastic_name'])
            self.annotation["description"] = self.annotation["mz"].apply(lambda x: clean_annotation[x]['re_name'])
            self.annotation["match_mz"] = self.annotation["mz"].apply(lambda x: clean_annotation[x]['match_mz'])
            print('Your spectra now has annotation')
        except:
            print('Problem with annotation, try to read or load spectra again')
    
    def save_spectra(self, name = 'spectra_name'):
        """Saves your spectra to database with chosen name"""
        con = sqlite3.connect("plastic_spectrum.db")
        cur = con.cursor()
        cur.execute('SELECT name FROM spectra WHERE name = ? LIMIT 1', (name,))
        if cur.fetchall() != []:
            cur.close()
            con.close()
            print('This name is already in database')
        else:
            data_to_load = self.annotation[["mz", "intensity", "residue_id", "status"]]
            data_to_load['name'] = name
            try:
                cur.executemany('INSERT INTO spectra(mz, intensity, residue_id, status, name) VALUES (?,?,?,?,?)', [tuple(i) for i in data_to_load.values])
                con.commit()
                cur.close()
                con.close()
                print('Successfully updated')
            except:
                cur.close()
                con.close()
                print('Update error')

    def load_spectra(self, name = 'spectra_name'):
        """Loads spectra from database based on name"""
        spec_name = name
        con = sqlite3.connect("plastic_spectrum.db")
        cur = con.cursor()
        try:
            query = cur.execute("SELECT mz, intensity, residue_id, status FROM spectra WHERE name = ?", (spec_name,))
            cols = [column[0] for column in query.description]
            loaded = pd.DataFrame.from_records(data = query.fetchall(), columns = cols)
            cur.close()
            con.close()
            loaded["plastic_name"] = loaded["residue_id"].apply( lambda x: self.__name_search(x)[0] if x>0 else None)
            loaded ["description"] = loaded["residue_id"].apply( lambda x: self.__name_search(x)[1] if x>0 else None)
            loaded ["match_mz"] = loaded["residue_id"].apply( lambda x: self.__name_search(x)[2] if x>0 else None)
            self.annotation = loaded
        except:
            print('Loading error')    
    
    def plot_spectra(self, title = 'Annotation results'):    
        """Plots your spectra"""
        exp_peaks = self.annotation.shape[0]
        theor_peaks = self.annotation[self.annotation['status'] != 'No matches'].shape[0]
        decoy_peaks = self.annotation[(self.annotation['status'] == 'Random match') | (self.annotation['status'] == 'Possibility of random match')].shape[0]
        match_peaks = self.annotation[self.annotation['status'] == 'Match'].shape[0]

        # experimental spectra
        plt.figure(figsize = (30,10))
        
        plt.bar(data = self.annotation,
                x = 'mz',
                height = 'intensity',
                color = "black",
                linewidth = 0.8)
        
        plt.xlabel(r'$m/z$', size = 15)
        plt.ylabel(r'$Intensity$', size = 15)
        plt.grid(True, axis = 'y', color = 'black', linestyle = ':', linewidth = 0.1)
        plt.title(f'{title}, experimental spectra: total experimental peaks — {exp_peaks}', size = 20)
        plt.show()
     
        # combined spectra random
        plt.figure(figsize = (30,10))
        
        plt.bar(data = self.annotation,
                x = 'mz',
                height = 'intensity',
                color = "black",
                linewidth = 0.8)
        plt.bar(data = self.annotation[self.annotation['status'].isin(['Match', 'Random match', 'Possibility of random match'])],
                x = 'match_mz',
                height =  'intensity',
                color = 'red',
                label = 'description',
                linewidth = 0.05)        
        
        plt.xlabel(r'$m/z$', size = 15)
        plt.ylabel(r'$Intensity$', size = 15)
        plt.grid(True, axis = 'y', color = 'black', linestyle = ':', linewidth = 0.1)
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.title(f'{title}, combined spectra with random matches: total experimental peaks — {exp_peaks}, total theoretical peaks — {theor_peaks}, total decoy peaks — {decoy_peaks}', size = 20)
        plt.show()

        # combined spectra nonrandom
        plt.figure(figsize = (30,10))
        
        plt.bar(data = self.annotation,
                x = 'mz',
                height = 'intensity',
                color = "black",
                linewidth = 0.8)
        plt.bar(data = self.annotation[self.annotation['status'] == 'Match'],
                x = 'match_mz',
                height = 'intensity',
                color = 'red',
                label = 'description',
                linewidth = 0.05)    
        
        plt.xlabel(r'$m/z$', size = 15)
        plt.ylabel(r'$Intensity$', size = 15)
        plt.grid(True, axis = 'y', color = 'black', linestyle = ':', linewidth = 0.1)
        plt.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.title(f'{title}, combined spectra without random matches: total experimental peaks — {exp_peaks}, total matching peaks — {match_peaks}', size = 20)
        plt.show()














class generated_plastic:


    def __init__(self, name):
        self.name = name
        self.residue = pd.DataFrame(columns = ["plastic_id", "name", "mz", "decoy_status"])
        # mono-isotopic masses
        self.elem_mass = {'C': 12.000000,
                        'H': 1.007825,
                        'S': 31.972071,
                        'O': 15.994915,
                        'Na': 22.989770,
                        'F': 18.998403,
                        'K': 39.0983}
        # electron mass
        self.e_mass = 0.00055
        # neutron mass
        self.n_mass = 1.008664
        # check and write plastic name if needed
        con = sqlite3.connect("plastic_spectrum.db")
        cur = con.cursor()
        query = cur.execute('SELECT plastic_id FROM plastic WHERE plastic_name = ?;', (name,))
        if query.fetchall() == []:
            cur.execute("INSERT INTO plastic(plastic_name) VALUES(?);", (name,))
            con.commit()
        else:
            print('plastic is already in database')
        cur.close()
        con.close()
        con = sqlite3.connect("plastic_spectrum.db")
        cur = con.cursor()
        query = cur.execute('SELECT plastic_id FROM plastic WHERE plastic_name = ?;', (name,))
        self.plastic_id = query.fetchall()[0][0]
        cur.close()
        con.close()
                


   def pss_gen(self, n = 60 , max_so3na = 5 , max_so3 = 5):
        """Generates PSS plastic database
            n - number of monomers
            max_so3na - number of SO3Na losses
            max_so3 -  number of SO3 losses"""
        self.name = self.name + '_' + str(n)
        elem_mass = self.elem_mass
        e_mass = self.e_mass
        n_mass = self.n_mass 
        
        # PSS mass calculation for actual polymer and its decoy
        polymer_mass = (elem_mass['C'] * 8 + elem_mass['H'] * 7 + elem_mass['S'] * 1 + elem_mass['O'] * 3 + elem_mass['Na'] * 1) * n 
        decoy_polym = polymer_mass + 1000

        # delta in case of 1% of isotops
        delta_mass = round(n_mass * 0.01 * 8 * n)
        
        charge_fragments = [] # list of masses for fragment ions
        fragment_labels = [] # list of information about fragment ions
        decoy_charge_fragments = [] # list of masses for decoy fragment ions
        decoy_fragment_labels = [] # list of information about decoy fragment ions
        small_fragments = []
        small_labels = []

        for i in range(1, n + 1): # list of possible na loss, depends on n
            for j in range(max_so3na + 1): # list of possible so3na loss
                for k in range(max_so3 + 1): # list of possible so3 loss
                    for l in range(21): # list of possible monomer units loss
                        for m in range(delta_mass - 2, delta_mass + 2): # list of possible isotopic error
                            if i + j > n or i + j - k <= 0: # sulphonate group lost but Na came back
                                continue
                            # mass loses calculations
                            na_loss = i * elem_mass['Na']
                            so3na_loss = j * (elem_mass['Na'] * 1 + elem_mass['S'] * 1 + elem_mass['O'] * 3)
                            so3_loss = k * (elem_mass['S'] * 1 + elem_mass['O'] * 3)
                            monom_loss = l * (elem_mass['C'] * 8 + elem_mass['H'] * 7 + elem_mass['S'] * 1 + elem_mass['O'] * 3 + elem_mass['Na'] * 1)
                            charge_state = i + j - k # charge state is equal to i + j - k
                            
                            # fragment ions calculations
                            fragment_mz = (polymer_mass + m - so3na_loss - na_loss + so3_loss - monom_loss + e_mass * charge_state)/charge_state
                            charge_fragments.append(fragment_mz) 
                            fragment_labels.append(f'PSS {n-l}-mer: {j} SO\u2083Na, {i} Na, {k} SO\u2083; Charge: {charge_state}; Peak(m/z): {fragment_mz}; Isotopic error: {m}')

                            # decoy fragment ions calculations
                            decoy_fragment_mz = (decoy_polym + m - so3na_loss - na_loss + so3_loss - monom_loss + e_mass * charge_state)/charge_state
                            decoy_charge_fragments.append(decoy_fragment_mz)
                            decoy_fragment_labels.append(f'PSS{n-l}-mer — decoy: {j} SO\u2083Na, {i} Na, {k} SO\u2083; Charge: {charge_state}; Peak(m/z): {fragment_mz}; Isotopic error: {m}')
                            
                            # small ions calculations
                            small_mz = - (m + so3na_loss + na_loss - so3_loss + monom_loss - e_mass * charge_state)/charge_state
                            small_fragments.append(small_mz)
                            small_labels.append(f'Charge: {j} SO\u2083Na, {i} Na, {k} SO\u2083; Charge: {charge_state}; Peak(m/z): {small_mz}; Isotopic error: {m}')

        labels = fragment_labels + decoy_fragment_labels + small_labels
        fragments = charge_fragments + decoy_charge_fragments + small_fragments

        generated_fragments = pd.DataFrame(zip(labels, fragments), columns = ['name', 'm/z'])
        generated_fragments['decoy_status'] =  generated_fragments['name'].apply(lambda x: 1 if 'decoy' in x else 0)
        self.residue[["name", "mz", "decoy_status"]] = generated_fragments[['name', 'm/z', 'decoy_status']]
        self.residue["plastic_id"] = self.residue["plastic_id"].fillna(self.plastic_id)



    def pfas_gen(self, n_c, n_f, n_k, n_s, n_o, n_h): 
        """Generates PFASes fragments database
        n_c - number of carbon
        n_f - number of fluoride
        n_k - number of potassium
        n_s - number of sulphur
        n_o - number of oxygen
        n_h - number of hydrogen"""
        # mono-isotopic masses
        elem_mass = self.elem_mass

        # electron mass
        e_mass = self.e_mass

        # polymer masses, depends on n
        pfas = elem_mass['C'] * n_c + elem_mass['H'] * n_h + elem_mass['S'] * n_s + elem_mass['O'] * n_o + elem_mass['K'] * n_k + elem_mass['F'] * n_f


        charge_fragments = [] # list of masses
        fragment_labels = [] # list of information about numbers of k (h) and SO3 (k) lost
        decoy_status = [] # decoy status


        if n_k != 0:
            for i in range(n_k + 1): # list of possible potassium loss
                for j in range(n_s + 1): # list of possible so3k loss
                    for k in range(n_s + 1): # list of possible so3 loss
                        if i + j - k <= 0: # sulphonate group lost but k came back
                            continue
                        k_loss = i * elem_mass['K']
                        so3k_loss = j * (elem_mass['K'] * 1 + elem_mass['S'] * 1 + elem_mass['O'] * 3)
                        so3_loss = k * (elem_mass['S'] * 1 + elem_mass['O'] * 3)
                        charge_state = i + j - k # charge state is equal to i + j - k
                        fragment_mz = (pfas - so3k_loss - k_loss + so3_loss + e_mass * charge_state)/charge_state
                        adduct_mass = pfas * 2 - i * elem_mass['K']
                        charge_fragments.append(fragment_mz)
                        charge_fragments.append(adduct_mass)
                        charge_fragments.append(so3k_loss) 
                        charge_fragments.append(so3_loss)
                        charge_fragments.append(k_loss)
                        fragment_labels.append(f'{i} K, {j} SO\u2083K, {k} SO\u2083 are lost; Charge: {charge_state}; Peak (m/z): {fragment_mz}')
                        decoy_status.append(0)
                        fragment_labels.append(f'Adduct\'s formula C{n_c * 2}F{n_f * 2}K{n_k}S{n_s * 2}O{n_o * 2}H{n_h * 2}: ; Peak (m/z): {adduct_mass}')
                        decoy_status.append(0)
                        fragment_labels.append(f'Fragment\'s formula K{n_k}S{n_s * 2}O{n_o * 2}:; Peak (m/z): {so3k_loss}')
                        decoy_status.append(0)
                        fragment_labels.append(f'Fragment\'s formula K{n_k}S{n_s * 2}O{n_o * 2}:; Peak (m/z): {so3_loss}')
                        decoy_status.append(0)
                        fragment_labels.append(f'Fragment\'s formula K{n_k}S{n_s * 2}O{n_o * 2}:; Peak (m/z): {k_loss}')
                        decoy_status.append(0)
        elif n_h != 0:
            for i in range(n_h + 1): # list of possible potassium loss
                if i <= 0: # sulphonate group lost but k came back
                    continue
                h_loss = i * elem_mass['H']
                charge_state = i
                fragment_mz = (pfas - h_loss + e_mass * charge_state)/charge_state
                adduct_mass = pfas * 2 - i * elem_mass['H']
                charge_fragments.append(fragment_mz)
                charge_fragments.append(adduct_mass)
                fragment_labels.append(f'{i} H are lost; Charge: {charge_state}; Peak (m/z): {fragment_mz}')
                decoy_status.append(0)
                fragment_labels.append(f'Adduct\'s formula C{n_c * 2}F{n_f * 2}O{n_o * 2}H{n_h * 2}; Charge: {charge_state}; Peak (m/z): {adduct_mass}')
                decoy_status.append(0)

        generated_fragments = pd.DataFrame(zip(fragment_labels, charge_fragments, decoy_status), columns = ['name', 'm/z', 'decoy_status'])
        generated_fragments = generated_fragments[generated_fragments['m/z'] != 0]
        self.residue[["name", "mz", "decoy_status"]] = generated_fragments[['name', 'm/z', 'decoy_status']]
        self.residue["plastic_id"] = self.residue["plastic_id"].fillna(self.plastic_id)
    
    
    def plastic_to_database(self):
        if len(self.residue.values) == 0:
            print('first_fill_residue')
        else:
            con = sqlite3.connect("plastic_spectrum.db")
            cur = con.cursor()
            try:
                cur.executemany('INSERT INTO residue(plastic_id, name, mz, decoy_status) VALUES(?, ?, ?, ?)', [tuple(i) for i in self.residue.values])
                print
            except:
                pass
            con.commit()
            cur.close()
            con.close()
        


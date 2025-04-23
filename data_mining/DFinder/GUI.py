import sys
sys.path.insert(0,'./modules')
from data_viewer import *
from data_reader import *
from data_writer import *
from duplicate_finder import *
from draw_utils import *
import PySimpleGUI as psg


def program():
	text = """INSTRUCTIONS
1. Chooce a csv file (clean *_log2.csv).
2. Select either masstime or value as criteria for duplicate search.
3. Set the appropriate precision.
4. After search, n duplicates are listed (0 to n-1).
5. Chooce duplicates by index for inspection.""".strip().split('\n')
	psg.popup(*text, line_width=70)
	# ask the caller to specify csv file
	fname = psg.popup_get_file('Select a csv file', title='File selector')
	try:
		if(fname.split('.')[1]!='csv'):
			psg.popup('Improper file, exiting...')
			return 1
	except:
		psg.popup('Improper file, exiting...')
		return 1
	
	# read the data
	data = read_file(fname)
	# chooce similarity measure
	layout0 = [[psg.Text('Criteria for similarity')],[psg.Button('masstime'), psg.Button('value')]]
	layout1 = [[psg.Text('Set similarity threshold\n(proportion of duplicate values')],
   [psg.Text('th: ', size=(5,1)), psg.Input(expand_x=True)],
   [psg.OK(), psg.Exit()]]
	layout2 = [
   [psg.Text('Set precision for mass and time similarity')],
   [psg.Text('dm: ', size=(5,1)),psg.Input(expand_x=True)],
   [psg.Text('dt: ', size=(5,1)), psg.Input(expand_x=True)],
   [psg.OK(), psg.Exit()]]
   
   
   
	window = psg.Window('Form', layout0)
	event, values = window.read()
	simtype = event
	if simtype=='masstime':
		window.close()
		window = psg.Window('Form', layout2)
	else:
		window.close()
		window = psg.Window('Form', layout1)
	
	event, values = window.read()
	if simtype=='masstime':
		indices, mt_pairs, similarities = get_MTime_duplicates(data, float(values[0]), float(values[1]))
	else:
		th = max(1.0/40, float(values[0]))
		indices, mt_pairs, similarities = get_value_duplicates(data, th)
	window.close()
	
	n=len(indices)
	write_duplicates(mt_pairs, similarities)
	f = open('tmp.txt')
	text = f.read()
	psg.popup_scrolled(text, title='Duplicates', font=("clean", 12),
	 size=(50,10), non_blocking=True, location=(400,500))
	
	while True:
		response = psg.popup_get_text('choose duplicates by integer (Duplicates list)')
		if response is None:
			break
		try:
			index = int(response)
			if(index>=n):
				psg.popup('Chooce values between {} and {}'.format(0,n-1))
				continue
		except:
			continue
		a, b = array_pair(data, indices[index])
		
		text = pair_string(a,b,indices[index],data)
		psg.popup_scrolled(text, title='Duplicate pair {}'.format(index), font=("Helvetica", 14),
		 size=(50,12), non_blocking=True, location=(1600,500))
		
		fig = get_figure(a,b)
		layout3 = [[psg.Text('Groups')],[psg.Canvas(key='-CANVAS-')],[psg.Button('Ok')]]
		window3 = psg.Window('Duplicate compounds', layout3, size=(500, 500), finalize=True, element_justification='center', font='Helvetica 18')
		tkcanvas = draw_figure(window3['-CANVAS-'].TKCanvas, fig)
		event, values = window3.read()
		window3.close()
	
	fname = psg.popup_get_text('Enter filename to save results or cancel')
	if fname is not None:
		if fname!='':
			fname = './results/' + fname
			make_copy(fname)
	
	cleanup()
	

def main():
	while not program():
		lo = [[psg.Text('Continue to select a new file or exit')],[psg.Button('Continue'), psg.Button('Exit')]]
		win = psg.Window('Continue/Exit',lo, finalize=True, modal=True, element_justification='c')
		event, values = win.read()
		if event == 'Exit':
			break
		win.close()


if __name__=='__main__':
	main()
	
	

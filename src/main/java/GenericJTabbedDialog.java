/*
 * TabbedPaneDemo.java requires one additional file:
 *   images/middle.gif.
 */

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.border.Border;

import ij.IJ;

public class GenericJTabbedDialog implements ActionListener {
	static final int   DIVIDE_COLUMNS_BY    = 2; // JTextField.setColumns() seem to make it twice wider?
	static final int   DEFAULT_STRING_WIDTH = 8;
    static final Color COMMENT_COLOR = new Color(0, 150, 0);
	static final int   COMMENT_PADY = 20; // make comment/message line wider

	//0: label:'LabelUI', input:'TextFieldUI', inp_units:'PanelUI'
	ArrayList<ArrayList<JComponent>> components = new ArrayList<ArrayList<JComponent>>(); // null component means a message, not am input
	ArrayList<ArrayList<JLabel>> labels = new ArrayList<ArrayList<JLabel>>();
	ArrayList<JComponent> tabs = new ArrayList<JComponent>();
	ArrayList<JButton> buttons = new ArrayList<JButton>();
	ArrayList<JScrollPane> scrollPanes = new ArrayList<JScrollPane>();
	private JFrame frame = new JFrame();
	public int width, height;
	private int read_tab=0,read_component=0; // current number of tab/component to be read next
	String result=null;
	private JDialog jd;
	public GenericJTabbedDialog(String title) {
		this(title, 600, 800);
	}
	public GenericJTabbedDialog(String title, int width, int height) {
//		final JFrame frame= new JFrame();
		jd = new JDialog(frame , title, true);
		components.add(new ArrayList<JComponent>());
		labels.add(new ArrayList<JLabel>());
		tabs.add (new JPanel(false)); // first (yet nameless tab)
		this.width = width;
		this.height = height;
	}

	public void addTab(String tab_name) {
		addTab(tab_name, "");

	}
	public void addTab(String tab_name, String tab_tooltip) {
		// See if there are any components for the current tab. If none (such as when none are yet added, just rename existing)
		// If there are some components already - start a new tab
		int ntab = tabs.size()-1;
		if (!components.get(ntab).isEmpty()) {
			components.add(new ArrayList<JComponent>());
			labels.add(new ArrayList<JLabel>());
			tabs.add (new JPanel(false)); // starting next tab
			ntab++;
		}
		tabs.get(ntab).putClientProperty("title",  tab_name);    // to be later used for tab name
		tabs.get(ntab).putClientProperty("tooltip",tab_tooltip); // to be later used for tab tooltip
	}

	public void addButtons(
			String [] labels,
			String [] actions,
			String [] tooltips) {

		for (int i = 0; i < labels.length; i++) {
			JButton button = new JButton(labels[i]);
			if ((actions != null) && (actions[i] != null)) button.setActionCommand(actions[i]);
			if ((tooltips != null) && (tooltips[i] != null)) button.setToolTipText(tooltips[i]);
			buttons.add(button);
		}
	}

	public void addDefaultButtons() {
		String [] labels = {"OK", "Cancel"};
		addButtons(labels, labels, labels);

	}

	// Common part of adding messages/input fields
	private void addLine(String prompt,
			             JComponent component, // may be null
			             String tooltip) {
		JLabel label;
   		if (prompt.indexOf('_') == -1) {
   			label = new JLabel(prompt);
   		} else {
   			label = new JLabel(prompt.replace('_', ' '));

   		}
		if (tooltip != null) {
			label.setToolTipText(tooltip);
			if (component != null) {
				component.setToolTipText(tooltip);
				Component [] comps = component.getComponents();
				if ((comps != null) && (comps.length >0)) {
					((JComponent) comps[0]).setToolTipText(tooltip);
				}
			}
		}
		labels.get(labels.size()-1).add(label);
		components.get(components.size()-1).add(component);
	}

	public void addMessage(String message) {
		addMessage(message, null);
	}
	public void addMessage(String message, String tooltip) {
		addLine(message, null, tooltip);
	}

	public void addStringField (String label, String value) {
		addStringField(label, value, DEFAULT_STRING_WIDTH, null);
	}


	public void addStringField (String label, String value,	String tooltip) {
		addStringField(label, value, DEFAULT_STRING_WIDTH, tooltip);
	}

	public void addStringField (String label, String value,	int width) {
		addStringField(label, value, width, null);
	}

	public void addStringField ( // no units here
			String label,
			String value,
			int    width,
			String tooltip) {
		JTextField inp = new JTextField();
		inp.setText(value);
		inp.setColumns(width/DIVIDE_COLUMNS_BY + 2); ///2);
		JPanel inp_units = new JPanel(false);
		inp_units.add(inp);
		inp_units.putClientProperty("type",  "String");
		inp_units.setLayout(new FlowLayout(FlowLayout.LEFT));
		addLine(label, inp_units, tooltip);

	}

	public void addNumericField(String label, double defaultValue, int digits) { // as in IJ
		addNumericField(label, defaultValue, digits, 6, null, null);
	}

	public void addNumericField(String label, double defaultValue, int digits, int columns, String units) { // as in IJ
		addNumericField(label, defaultValue, digits, columns, units, null);
	}

	public void addNumericField(String label, double defaultValue, int digits, int columns, String units, String tooltip) {
		String svalue = IJ.d2s(defaultValue, digits); // , columns
		JTextField inp = new JTextField();
		inp.setText(svalue);
		inp.setColumns(columns/DIVIDE_COLUMNS_BY + 2);
		JPanel inp_units = new JPanel(false);
		inp_units.add(inp);
		if ((units != null) && !units.isEmpty()) {
			JLabel junits = new JLabel(units);
			Font fnt = junits.getFont();
			junits.setFont(fnt.deriveFont(fnt.getStyle() | Font.ITALIC));
			inp_units.add(junits);
		}
		inp_units.putClientProperty("type",  "double");
		inp_units.setLayout(new FlowLayout(FlowLayout.LEFT));
		addLine(label, inp_units, tooltip);
//		Component[] comps = inp_units.getComponents();
	}

	public void addCheckbox(String label, boolean defaultValue) {
		addCheckbox(label, defaultValue, null);
	}

    public void addCheckbox(String label, boolean defaultValue, String tooltip) { // no support for preview functionality
    	JCheckBox checkbox = new JCheckBox(null, null, defaultValue); // text, icon, selected
		checkbox.putClientProperty("type",  "boolean");
    	addLine(label, checkbox, tooltip);
    }
 // Add choice too.
    public void addChoice(String label, String[] items, String defaultItem) {
    	addChoice(label, items, defaultItem, null);
    }
    public void addChoice(String label, String[] items, String defaultItem, String tooltip) {
    	int index = 0;
    	if (defaultItem != null) {
    		for (int i = 0; i < items.length; i++) if (items[i].equals(defaultItem)) {
    			index = i;
    			break;
    		}
    	}
    	JComboBox<String> combo = new JComboBox<String>(items);
    	combo.setSelectedIndex(index);

		JPanel inp_units = new JPanel(false);
		inp_units.add(combo);

		inp_units.putClientProperty("type",  "combo");
//		combo.setPreferredSize(new Dimension(200));
//		combo.setSize(200, combo.getPreferredSize().height);
//		combo.setMaximumSize(20); //  combo.getPreferredSize() );
		inp_units.setLayout(new FlowLayout(FlowLayout.LEFT));

    	addLine(label, inp_units, tooltip);
    }


    public void buildDialog() { // non-blocking, does not show

    	if (buttons.isEmpty()) {
    		addDefaultButtons();
    	}

		jd.setLayout(new BorderLayout()); //  new GridLayout(2, 2) );
		// prepare each panel first, then optionally make a tabbed layout
		for (int ntab = 0; ntab < tabs.size(); ntab++) {
			JComponent            tab = tabs.get(ntab);
			ArrayList<JComponent> tabComponents = components.get(ntab);
			ArrayList<JLabel>     tabLabels = labels.get(ntab);
	        tab.setLayout(new GridBagLayout());
	        GridBagConstraints gbc = new GridBagConstraints();
	        gbc.fill = GridBagConstraints.HORIZONTAL;
	    	gbc.weightx = 0.5;
	    	gbc.weighty = 0.0;
			for (int ncomp = 0; ncomp < tabLabels.size(); ncomp++) {
				JLabel label =       tabLabels.get(ncomp);
				JComponent component = tabComponents.get(ncomp);
	        	gbc.gridx = 0;
	        	gbc.gridy = ncomp;
	        	gbc.gridwidth = 1;
	        	gbc.ipady = 0;
				if (component == null) { // message/comment line
					label.setHorizontalAlignment(JLabel.CENTER);
					label.setForeground(COMMENT_COLOR);
	        		Font fnt = label.getFont();
	        		label.setFont(fnt.deriveFont(fnt.getStyle() | Font.ITALIC));
	        		gbc.gridwidth = 2;
	            	gbc.ipady = COMMENT_PADY;
	            	tab.add(label,gbc);
				} else {
					label.setHorizontalAlignment(JLabel.RIGHT);
	            	tab.add(label,gbc);
	        		gbc.gridx = 1;
	        		tab.add(component,gbc);
				}
			}
			// Add empty label that would push up all other rows
	        JLabel label_push = new JLabel("");
	    	gbc.weightx = 0.5;
	    	gbc.weighty = 0.5;
	    	gbc.gridx = 0;
	    	gbc.gridy = tabLabels.size();
	    	tab.add(label_push,gbc);

	    	// Wrap each panel with scroll panel
	        JScrollPane scrollPane = new JScrollPane(tab);
	        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED); // HORIZONTAL_SCROLLBAR_NEVER);
	        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
	    	scrollPanes.add(scrollPane);
		}
        Border padding = BorderFactory.createEmptyBorder(20,20,5,20);
		if (tabs.size() == 1) {
	        jd.add(new JScrollPane(scrollPanes.get(0)), BorderLayout.CENTER); // no tabs
		} else { // tabbed
			JTabbedPane tabbedPane = new JTabbedPane();
			for (int ntab = 0; ntab < tabs.size(); ntab++) {
				JScrollPane scrollPane = scrollPanes.get(ntab);
				JComponent tab = tabs.get(ntab);
				scrollPane.setBorder(padding);
				tabbedPane.addTab(
						(String) tab.getClientProperty("title"),
						null, // icon
						scrollPane,
						(String) tab.getClientProperty("tooltip"));
			}
	        jd.add(tabbedPane, BorderLayout.CENTER); //
		}
		// Add buttons
        JPanel panelButtons = new JPanel(false);
        panelButtons.setLayout(new FlowLayout());
        for (JButton b : buttons) {
        	panelButtons.add(b);
        	 b.addActionListener(this);
        }
        jd.add(panelButtons, BorderLayout.PAGE_END);
        jd.setSize(width,height);
    }

        public String showDialogAny() {
  			buildDialog();
        	jd.setVisible(true);
            return result;
        }
        public boolean showDialog() {
  			buildDialog();
        	jd.setVisible(true);
        	if ((result != null) && !result.equals("Cancel")) return true;
            return false;
        }

    @Override
	public void actionPerformed(ActionEvent e) {
    	System.out.println(e.getActionCommand());
    	if (e.getActionCommand().equals("Cancel")) {
    		result = e.getActionCommand();
    		jd.dispose();
    	} else if (e.getActionCommand().equals("OK")) {
    		result = e.getActionCommand();
    		jd.dispose();
    	}

    }

    public boolean wasCanceled() {
    	return (result == null) || (result.equals("Cancel"));
    }

    private boolean skipComments (
    		boolean incTab) { // increment tab if nothing left in this one (false - keep read_component equal to tab's components length
    	boolean got_it = false;
//    	boolean nothing_left = false;
//    	read_component ++;
    	while (!got_it) {
    		if (read_component >= labels.get(read_tab).size()) {
    			if (!incTab || (read_tab >= labels.size())) break;
    			read_tab++;
    			read_component = 0;
    		}
    		if (components.get(read_tab).get(read_component) != null) {
    			got_it = true;
    		} else {
    			read_component++;
    		}
    	}
    	return got_it;
    }

    public void checkType (String type) { // type == null for null components OK
    	JComponent component = components.get(read_tab).get(read_component);
    	String this_type = "null";
    	if (component != null) {
    		this_type = (String) component.getClientProperty("type");
    		if (this_type.equals(type)) return;
    	} else if (type == null) return;
    	// component type mismatch - report
    	String msg="Dialog components mismatch for tab="+read_tab+", component="+read_component+
    			", label='"+labels.get(read_tab).get(read_component).getText()+"': expected "+type+", got "+this_type;
    	IJ.showMessage(msg);
    	throw new IllegalArgumentException (msg);
    }

    public boolean getNextBoolean() {
    	skipComments(true);
    	checkType("boolean");
    	boolean rslt = ((JCheckBox) components.get(read_tab).get(read_component)).isSelected();
    	read_component ++;
    	return rslt;

    }

    public String getNextString() {
    	skipComments(true); // add testing for checkbox?
    	checkType("String");
    	String rslt =  ((JTextField) components.get(read_tab).get(read_component).getComponents()[0]).getText();
    	read_component ++;
    	return rslt;
    }

    public double getNextNumber() {
    	skipComments(true); // add testing for checkbox?
    	checkType("double");
    	double rslt =  Double.parseDouble(((JTextField) components.get(read_tab).get(read_component).getComponents()[0]).getText());
    	read_component ++;
    	return rslt;
    }

    public int getNextChoiceIndex() {
    	skipComments(true); // add testing for checkbox?
    	checkType("combo");
    	@SuppressWarnings("unchecked")
		JComboBox<String> combo = (JComboBox<String>) components.get(read_tab).get(read_component).getComponents()[0];
    	String selectedItem = (String)combo.getSelectedItem();
//    	Object [] items = combo.getSelectedObjects();
    	int size = combo.getItemCount();
    	int index = 0;
    	if (selectedItem != null) {
    		for (int i = 0; i < size; i++) if (combo.getItemAt(i).equals(selectedItem)) {
    			index = i;
    			break;
    		}
    	}


    	System.out.println(combo.getSelectedItem());
//    	for (int i = 0; i < selectedItem.length; i++){
 //   	    System.out.println(String.format("item %s = %s", i, selectedItem[i]));
  //  	}

    	read_component ++;
    	return index;
    }

/*
String[] selectedItem = (String[])combo.getSelectedItem();
for (int i = 0; i < selectedItem.length; i++){
    System.out.println(String.format("item %s = %s", i, selectedItem[i]));
}
 */

}
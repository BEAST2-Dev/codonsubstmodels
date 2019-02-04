package beast.app.beauti;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.GeneticCode;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Map;


public class CodonAlignmentViewer extends JPanel {
    private static final long serialVersionUID = 1L;

    Object[][] tableData;
    Object[] columnData;
    boolean useColor = false;
    // flag to indicate that the most frequently occurring character is shown as a dot
    boolean useDots = true;
    CodonAlignment m_alignment;
    Map<Character, Color> m_customColorMap = new HashMap<>();

    JTable mainTable;

    private int stopCodonSite = -1;

    /**
     * constructor processes alignment and sets up table with first column fixed *
     */
    public CodonAlignmentViewer(Alignment data) {
        if (! (data instanceof CodonAlignment) )
            throw new IllegalArgumentException("Codon alignment is required ! ");
        m_alignment = (CodonAlignment) data;

        int siteCount = data.getSiteCount();
        int taxonCount = data.getTaxonCount();
        tableData = new Object[taxonCount][siteCount + 1];
        char[] headerChar = updateTableData();

        // set up row labels
        for (int i = 0; i < taxonCount; i++) {
            tableData[i][0] = data.getTaxaNames().get(i);
        }

        // set up column labels
        columnData = new Object[siteCount + 1];
        updateColumnData(headerChar);

        // create table in scrollpane with first column fixed
        final TableModel fixedColumnModel = new AbstractTableModel() {
            private static final long serialVersionUID = 1L;

            @Override
            public int getColumnCount() {
                return 1;
            }

            @Override
            public String getColumnName(int column) {
                return columnData[column] + "";
            }

            @Override
            public int getRowCount() {
                return tableData.length;
            }

            @Override
            public Object getValueAt(int row, int column) {
                return tableData[row][column];
            }
        };

        final TableModel mainModel = new AbstractTableModel() {
            private static final long serialVersionUID = 1L;

            @Override
            public int getColumnCount() {
                return columnData.length - 1;
            }

            @Override
            public String getColumnName(int column) {
                return columnData[column + 1] + "";
            }

            @Override
            public int getRowCount() {
                return tableData.length;
            }

            @Override
            public Object getValueAt(int row, int column) {
                return tableData[row][column + 1];
            }
        };

        JTable fixedTable = new JTable(fixedColumnModel);
        fixedTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        Font font = fixedTable.getFont();
        font = new Font(font.getFontName(), font.getStyle(), font.getSize() * 2/3);
        fixedTable.setFont(font);
        TableColumn col = fixedTable.getColumnModel().getColumn(0);
        col.setPreferredWidth(200);
        fixedTable.getTableHeader().setFont(font);

        mainTable = new JTable(mainModel);
        mainTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        mainTable.setFont(font);
        mainTable.getTableHeader().setFont(font);
        for (int i = 0; i < siteCount; i++) {
            col = mainTable.getColumnModel().getColumn(i);
            col.setPreferredWidth(6);
        }

        ListSelectionModel model = fixedTable.getSelectionModel();
        mainTable.setSelectionModel(model);
        mainTable.setShowGrid(false);

        JScrollPane scrollPane = new JScrollPane();
        scrollPane.getViewport().add(mainTable);
        Dimension fixedSize = fixedTable.getPreferredSize();
        JViewport viewport = new JViewport();
        viewport.setView(fixedTable);
        viewport.setPreferredSize(fixedSize);
        viewport.setMaximumSize(fixedSize);
        scrollPane.setCorner(ScrollPaneConstants.UPPER_LEFT_CORNER, fixedTable.getTableHeader());
        scrollPane.setRowHeaderView(viewport);


        setLayout(new BorderLayout());
        add(scrollPane, BorderLayout.CENTER);
    }

    private void updateColumnData(char[] headerChar) {
        // set up column labels
        for (int i = 1; i < columnData.length; i++) {
            columnData[i] = "<html>.<br>" + headerChar[i - 1] + "</html>";
        }
        columnData[0] = "<html><br>taxon name</html>";
        columnData[1] = "<html>1<br>" + headerChar[0] + "</html>";
        for (int i = 10; i < columnData.length-1; i += 10) {
            String s = i + "";
            for (int j = 0; j < s.length(); j++) {
                if (i+j < columnData.length) {
                    columnData[i + j] = "<html>" + s.charAt(j) + "<br>" + headerChar[i - 1] + "</html>";
                }
            }
            columnData[i - 5] = "<html>+<br>" + headerChar[i - 1] + "</html>";
        }

    }

    private char[] updateTableData() {
        int siteCount = m_alignment.getSiteCount();
        int taxonCount = m_alignment.getTaxonCount();

        // set up table content
        DataType dataType = m_alignment.getDataType();
//        if (! (dataType instanceof Codon) )
//            throw new IllegalArgumentException("Codon data type is required ! " + dataType.getTypeDescription());

        char[] headerChar = new char[siteCount];
        Object[][] colorMap = setupColorMap();

        stopCodonSite = -1;
        try {
            for (int i = 0; i < siteCount; i++) {
                int patternIndex_ = m_alignment.getPatternIndex(i);
                int[] pattern = m_alignment.getPattern(patternIndex_);
//                String patternString = dataType.state2string(pattern);
                String patternString = ((Codon)dataType).stateToAminoAcid(pattern);

                headerChar[i] = mostFrequentCharInPattern(patternString);
                for (int j = 0; j < taxonCount; j++) {
                    char c = patternString.charAt(j);
                    if (c == headerChar[i]) {
                        tableData[j][i + 1] = colorMap[0][c];
                    } else {
                        tableData[j][i + 1] = colorMap[1][c];
                    }
                }

                if (stopCodonSite < 0 && patternString.contains(Character.toString(Codon.STOP_CHARACTER))) {
                    stopCodonSite = i;
                }
            }
        } catch (Exception e) {
            // ignore
        }
        return headerChar;
    }

    /**
     * determine content of table cells.
     * Without color, only Characters are displayed, which can be a bit faster than using color
     * With color, the color is encoded in HTML
     *
     * @return an array of 2x256 where the first entry is for the most frequently occurring character,
     *         and the second for the others
     *         *
     */
    private Object[][] setupColorMap() {
        if (useColor) {
            String[][] colorMap = new String[2][256];
            for (int k = 'A'; k < 'Z'; k++) {
                int i = k - 'A';
                int red = ((i & 0x80) >> 7) + ((i & 0x10) >> 4) + ((i & 0x2) << 1);
                int green = ((i & 0x40) >> 6) + ((i & 0x08) >> 2) + ((i & 0x4));
                int blue = ((i & 0x20) >> 5) + ((i & 0x04) >> 1) + ((i & 0x1) << 2);
                int color = (red << 21 + (green << 18)) + (green << 13) + (blue << 10) + (blue << 5) + (red << 2);
                colorMap[0][k] = "<html><font color='#" + Integer.toString(color, 16) + "'><b>.</b></html>";
                colorMap[1][k] = "<html><font color='#" + Integer.toString(color, 16) + "'><b>" + ((char) k) + "</font></html>";
            }
            for (char c : m_customColorMap.keySet()) {
                Color color = m_customColorMap.get(c);
                colorMap[0][c] = "<html><font color='#" + Integer.toString(color.getRGB(), 16) + "'><b>.</b></html>";
                colorMap[1][c] = "<html><font color='#" + Integer.toString(color.getRGB(), 16) + "'><b>" + c + "</font></html>";
            }
            if (!this.useDots) {
                colorMap[0] = colorMap[1];
            }
            return colorMap;
        } else {
            Character[][] colorMap = new Character[2][256];
            for (int i = 0; i < 256; i++) {
                colorMap[0][i] = '.';
                colorMap[1][i] = (char) i;
            }
            if (!this.useDots) {
                colorMap[0] = colorMap[1];
            }
            return colorMap;
        }
    }

    private char mostFrequentCharInPattern(String pattern) {
        char[] counts = new char[256];
        for (int i = 0; i < pattern.length(); i++) {
            counts[pattern.charAt(i)]++;
        }
        int maxIndex = 0, max = 0;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > max) {
                maxIndex = i;
                max = counts[i];
            }
        }
        return (char) maxIndex;
    }

    public void showInDialog() {
        JDialog dlg = new JDialog();
        dlg.setName("CodonAlignmentViewer");
        dlg.add(this);

        Box buttonBox = Box.createHorizontalBox();
        JCheckBox useDotsCheckBox = new JCheckBox("Use dots", true);
        useDotsCheckBox.addActionListener(e -> {
            JCheckBox _useDots = (JCheckBox) e.getSource();
            useDots = _useDots.isSelected();
            updateTableData();
            repaint();
        });
        buttonBox.add(useDotsCheckBox);

        JCheckBox useColorCheckBox = new JCheckBox("Use Color");
        useColorCheckBox.setName("UseColor");
        useColorCheckBox.addActionListener(e -> {
            JCheckBox hasColor = (JCheckBox) e.getSource();
            useColor = hasColor.isSelected();
            updateTableData();
            repaint();
        });
        buttonBox.add(useColorCheckBox);

        JComboBox geneticCodeComboBox = new JComboBox(GeneticCode.GENETIC_CODE_NAMES);
        int idx = 0;
        for (int i = 0; i < GeneticCode.GENETIC_CODE_NAMES.length; i++) {
            String gcName = GeneticCode.GENETIC_CODE_NAMES[i];
            if (gcName.equalsIgnoreCase(m_alignment.getGeneticCode().getName())) {
                idx = i;
                break;
            }
        }
        geneticCodeComboBox.setSelectedIndex(idx);
        geneticCodeComboBox.addActionListener(e -> {
            JComboBox cb = (JComboBox)e.getSource();
            String geneticCodeName = (String)cb.getSelectedItem();
            GeneticCode geneticCode = GeneticCode.findByName(geneticCodeName);
            m_alignment.setGeneticCode(geneticCode);

            char[] headerChar = updateTableData();
//            System.out.println(headerChar);
            updateColumnData(headerChar); // not refresh
            JTableHeader th = mainTable.getTableHeader();
            TableColumnModel tcm = th.getColumnModel();
            for (int i = 0; i < tcm.getColumnCount(); i++) {
                TableColumn tc = tcm.getColumn(i);
                tc.setHeaderValue( columnData[i+1] );
            }
            th.repaint();
            repaint();
        });
        buttonBox.add(geneticCodeComboBox);

        JButton checkStopCodonJButton = new JButton("Check Stop Codon");
        checkStopCodonJButton.addActionListener(e -> {
            if (stopCodonSite < 0) {
                JOptionPane.showMessageDialog(this,
                        "There is no stop codon,\nthe genetic code is valid for this alignment.",
                        "No stop codon", JOptionPane.INFORMATION_MESSAGE);
            } else {
                JOptionPane.showMessageDialog(this,
                        "Find a stop codon at triplet " + (stopCodonSite+1) +
                                ",\nplease choose the correct genetic code or use codon alignment!",
                        "Find stop codon", JOptionPane.ERROR_MESSAGE);
            }
        });
        buttonBox.add(checkStopCodonJButton);

        dlg.add(buttonBox, BorderLayout.SOUTH);

        int size = UIManager.getFont("Label.font").getSize();
        dlg.setSize(1024 * size / 13, 600 * size / 13);
        dlg.setModal(true);
        dlg.setVisible(true);
        dlg.dispose();
    }

}

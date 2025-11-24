#!/usr/bin/env python3
"""
GSL + Ceramide Transition Generator - GUI Application
A graphical interface for generating lipid transition lists
"""

__author__ = "Andreas J. Hülsmeier"
__copyright__ = "Copyright 2025, Andreas J. Hülsmeier / University of Zurich, University Hospital Zurich"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Andreas J. Hülsmeier"
__email__ = "andreas.huelsmeier@uzh.ch"
__status__ = "Prototype"

import sys
import pandas as pd
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QComboBox,
    QCheckBox, QHBoxLayout, QPushButton, QTableWidget,
    QTableWidgetItem, QFileDialog, QMessageBox, QGroupBox,
    QLineEdit, QScrollArea, QTextEdit, QDialog, QListWidget,
    QAbstractItemView, QSpinBox
)
from PySide6.QtGui import QIcon
import gslgen_v10
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)

class ConfigDialog(QDialog):
    """Dialog for configuring LCB and fatty acid selections"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Configure Lipid Chains")
        self.setMinimumSize(600, 500)

        config = gslgen_v10.ConfigManager.load_config()

        layout = QVBoxLayout()

        # ========== LCB SELECTION ==========

                # LCB Selection for Standard Ceramides
        layout.addWidget(QLabel("Standard Ceramide LCBs (2 OH) - Select multiple:"))
        self.lcb_list_standard = QListWidget()
        self.lcb_list_standard.setSelectionMode(QAbstractItemView.MultiSelection)

        standard_lcbs = ["16:0;2", "16:1;2", "17:0;2", "17:1;2",
                         "18:0;2", "18:0;3", "18:1;2", "18:2;2",
                         "19:0;2", "19:1;2", "20:0;2", "20:1;2"]

        for lcb in standard_lcbs:
            self.lcb_list_standard.addItem(lcb)

        for i in range(self.lcb_list_standard.count()):
            item = self.lcb_list_standard.item(i)
            if item.text() in config["lcb_selections"]["standard"]:
                item.setSelected(True)

        layout.addWidget(self.lcb_list_standard)

        # LCB Selection for Deoxy Ceramides
        layout.addWidget(QLabel("\n1-Deoxy Ceramide LCBs (1 OH) - Select multiple:"))
        self.lcb_list_doxcer = QListWidget()
        self.lcb_list_doxcer.setSelectionMode(QAbstractItemView.MultiSelection)

        doxcer_lcbs = ["16:0;1", "16:1;1", "17:0;1", "17:1;1",
                       "18:0;1", "18:1;1", "18:2;1",
                       "19:0;1", "19:1;1", "20:0;1", "20:1;1"]

        for lcb in doxcer_lcbs:
            self.lcb_list_doxcer.addItem(lcb)

        for i in range(self.lcb_list_doxcer.count()):
            item = self.lcb_list_doxcer.item(i)
            if item.text() in config["lcb_selections"]["doxCer"]:
                item.setSelected(True)

        layout.addWidget(self.lcb_list_doxcer)


        # ========== FATTY ACID RANGE ==========
        layout.addWidget(QLabel("\nFatty Acid Chain Length:"))

        fa_range_layout = QHBoxLayout()
        fa_range_layout.addWidget(QLabel("Min:"))
        self.fa_min_spin = QSpinBox()
        self.fa_min_spin.setRange(12, 30)
        self.fa_min_spin.setValue(config["fatty_acid_range"]["min_length"])
        fa_range_layout.addWidget(self.fa_min_spin)

        fa_range_layout.addWidget(QLabel("Max:"))
        self.fa_max_spin = QSpinBox()
        self.fa_max_spin.setRange(12, 30)
        self.fa_max_spin.setValue(config["fatty_acid_range"]["max_length"])
        fa_range_layout.addWidget(self.fa_max_spin)

        layout.addLayout(fa_range_layout)

        # Even-chain only checkbox
        self.even_chain_checkbox = QCheckBox("Even-chain length only (e.g., 16, 18, 20...)")
        if config["fatty_acid_range"].get("even_chain_only", False):
            self.even_chain_checkbox.setChecked(True)
        layout.addWidget(self.even_chain_checkbox)

        # ========== FATTY ACID UNSATURATIONS (ONLY ONE INSTANCE) ==========
        layout.addWidget(QLabel("\nFatty Acid Unsaturations:"))
        self.unsat_checkboxes = {}
        unsat_layout = QHBoxLayout()
        for unsat in range(0, 4):
            checkbox = QCheckBox(f"{unsat} double bond{'s' if unsat != 1 else ''}")
            if unsat in config["fatty_acid_range"]["unsaturations"]:
                checkbox.setChecked(True)
            self.unsat_checkboxes[unsat] = checkbox
            unsat_layout.addWidget(checkbox)
        layout.addLayout(unsat_layout)

        # Add visual separator
        separator = QLabel()
        separator.setStyleSheet("border-top: 1px solid #ccc; margin: 10px 0;")
        layout.addWidget(separator)

        # Create preview label with styling
        self.preview_label = QLabel("Expected combinations: Calculating...")
        self.preview_label.setStyleSheet(
            "color: #0066cc; "
            "font-weight: bold; "
            "font-size: 11pt; "
            "padding: 8px; "
            "background-color: #f0f8ff; "
            "border-radius: 4px;"
        )
        layout.addWidget(self.preview_label)

        # Connect all widgets to update the preview
        self.lcb_list_standard.itemSelectionChanged.connect(self.update_preview)
        self.lcb_list_doxcer.itemSelectionChanged.connect(self.update_preview)
        self.fa_min_spin.valueChanged.connect(self.update_preview)
        self.fa_max_spin.valueChanged.connect(self.update_preview)
        self.even_chain_checkbox.toggled.connect(self.update_preview)
        for checkbox in self.unsat_checkboxes.values():
            checkbox.toggled.connect(self.update_preview)

        # Calculate initial preview
        self.update_preview()

        # ========== INFO AND BUTTONS ==========
        info_label = QLabel("\nNote: Changes will be saved and used for all future generations.")
        info_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(info_label)

        button_layout = QHBoxLayout()
        save_btn = QPushButton("Save Configuration")
        cancel_btn = QPushButton("Cancel")
        save_btn.clicked.connect(self.save_and_close)
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(save_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def save_and_close(self):
        """Save configuration and close dialog"""
        min_length = self.fa_min_spin.value()
        max_length = self.fa_max_spin.value()

        if max_length < min_length:
            QMessageBox.warning(
                self,
                "Invalid Range",
                f"Maximum fatty acid length ({max_length}) must be >= "
                f"minimum length ({min_length})!"
            )
            return

        # Get selected LCBs from both lists
        selected_standard = [self.lcb_list_standard.item(i).text()
                            for i in range(self.lcb_list_standard.count())
                            if self.lcb_list_standard.item(i).isSelected()]

        selected_doxcer = [self.lcb_list_doxcer.item(i).text()
                          for i in range(self.lcb_list_doxcer.count())
                          if self.lcb_list_doxcer.item(i).isSelected()]

        # Check if at least one LCB is selected
        if not selected_standard and not selected_doxcer:
            QMessageBox.warning(self, "No Selection", "Please select at least one LCB!")
            return

        # Get unsaturations
        selected_unsaturations = [unsat for unsat, checkbox in self.unsat_checkboxes.items()
                                 if checkbox.isChecked()]

        if not selected_unsaturations:
            QMessageBox.warning(self, "No Selection", "Please select at least one unsaturation level!")
            return

        # Build config with separate standard and doxCer lists
        config = {
            "lcb_selections": {
                "standard": selected_standard,
                "doxCer": selected_doxcer
            },
            "fatty_acid_range": {
                "min_length": self.fa_min_spin.value(),
                "max_length": self.fa_max_spin.value(),
                "unsaturations": selected_unsaturations,
                "even_chain_only": self.even_chain_checkbox.isChecked()
            },
            "selected_fatty_acids": None
        }

        # Save
        gslgen_v10.ConfigManager.save_config(config)
        QMessageBox.information(self, "Saved", "Configuration saved successfully!")
        self.accept()

    def update_preview(self):
        """Calculate and display expected number of combinations"""
        try:
            # Count selected LCBs
            lcb_standard_count = sum(
                1 for i in range(self.lcb_list_standard.count())
                if self.lcb_list_standard.item(i).isSelected()
            )
            lcb_doxcer_count = sum(
                1 for i in range(self.lcb_list_doxcer.count())
                if self.lcb_list_doxcer.item(i).isSelected()
            )
            total_lcb_count = lcb_standard_count + lcb_doxcer_count

            # Calculate fatty acid count
            fa_min = self.fa_min_spin.value()
            fa_max = self.fa_max_spin.value()

            if fa_max < fa_min:
                fa_count = 0
                range_text = "(invalid range)"
            elif self.even_chain_checkbox.isChecked():
                fa_count = len(range(fa_min, fa_max + 1, 2))
                range_text = f"(even only: {fa_min}-{fa_max})"
            else:
                fa_count = fa_max - fa_min + 1
                range_text = f"({fa_min}-{fa_max})"

            # Count selected unsaturations
            unsat_count = sum(
                1 for checkbox in self.unsat_checkboxes.values()
                if checkbox.isChecked()
            )

            # Calculate total
            total = total_lcb_count * fa_count * unsat_count

            # Display with appropriate styling
            if total == 0:
                self.preview_label.setText(
                    "⚠️ No combinations will be generated - check your selections"
                )
                self.preview_label.setStyleSheet(
                    "color: #cc6600; font-weight: bold; font-size: 11pt; "
                    "padding: 8px; background-color: #fff8e6; border-radius: 4px;"
                )
            else:
                self.preview_label.setText(
                    f"📊 Expected combinations: {total:,}\n"
                    f"({total_lcb_count} LCBs × {fa_count} FAs {range_text} × "
                    f"{unsat_count} unsaturations)"
                )
                self.preview_label.setStyleSheet(
                    "color: #0066cc; font-weight: bold; font-size: 11pt; "
                    "padding: 8px; background-color: #f0f8ff; border-radius: 4px;"
                )
        except Exception:
            self.preview_label.setText("Preview unavailable")



class GSLGui(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GSL & Ceramide Transition Generator v1.0")
        self.setMinimumSize(1200, 900)

        # Initialize charge/adduct tracking dictionaries
        self.charge_checkboxes = {}
        self.adduct_checkboxes_by_charge = {}
        self.adduct_group_boxes = {}

        # Initialize auto-save blocking flag
        self.block_auto_save = False

        # Set window icon
        self.set_window_icon()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        main_widget = QWidget()
        layout = QVBoxLayout(main_widget)

        # Lipid Class Selection
        self.lipid_label = QLabel("Lipid Class:")
        self.lipid_combo = QComboBox()
        self.lipid_combo.addItems(gslgen_v10.LipidDatabase.get_all_classes())
        self.lipid_combo.currentTextChanged.connect(self.update_info)
        layout.addWidget(self.lipid_label)
        layout.addWidget(self.lipid_combo)

        # Info display
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(80)
        layout.addWidget(QLabel("Lipid Information:"))
        layout.addWidget(self.info_text)

        # === CHARGE STATES AND ADDUCT SELECTION ===
        selection_group = QGroupBox("Charge States and Adducts")
        selection_layout = QVBoxLayout()

        # Charge state selection
        charge_label = QLabel("Select Charge States:")
        charge_label.setStyleSheet("font-weight: bold;")
        selection_layout.addWidget(charge_label)

        charge_layout = QHBoxLayout()
        for charge in [1, 2, 3, 4, 5]:
            cb = QCheckBox(f"Charge ±{charge}")
            cb.setChecked(charge == 1)
            cb.stateChanged.connect(self.update_adduct_visibility)
            self.charge_checkboxes[charge] = cb
            charge_layout.addWidget(cb)
        charge_layout.addStretch()

        self.recommended_btn = QPushButton("Use Recommended")
        self.recommended_btn.setToolTip("Set recommended charge states for selected lipid class")
        self.recommended_btn.clicked.connect(self.set_recommended_charges)
        charge_layout.addWidget(self.recommended_btn)

        selection_layout.addLayout(charge_layout)

        # Adduct selection header
        adduct_label = QLabel("Select Adducts (for selected charge states):")
        adduct_label.setStyleSheet("font-weight: bold; margin-top: 10px;")
        selection_layout.addWidget(adduct_label)

        # Define adducts by charge state
        adducts_by_charge = {
            1: {'positive': ['[M+H]+', '[M+Na]+', '[M+NH4]+'],
                'negative': ['[M-H]-', '[M+CH3COO]-', '[M+HCOO]-']},
            2: {'positive': ['[M+2H]2+', '[M+H+Na]2+', '[M+2Na]2+'],
                'negative': ['[M-2H]2-']},
            3: {'positive': ['[M+3H]3+', '[M+2H+Na]3+', '[M+H+2Na]3+', '[M+3Na]3+'],
                'negative': ['[M-3H]3-']},
            4: {'positive': ['[M+4H]4+'],
                'negative': ['[M-4H]4-']},
            5: {'positive': ['[M+5H]5+'],
                'negative': ['[M-5H]5-']}
        }

        # Create adduct groups for each charge state
        for charge in [1, 2, 3, 4, 5]:
            group_box = QGroupBox(f"Charge ±{charge}")
            group_layout = QVBoxLayout()
            self.adduct_checkboxes_by_charge[charge] = {}

            # Positive adducts
            pos_label = QLabel("Positive Mode:")
            pos_label.setStyleSheet("font-style: italic; color: #0066cc;")
            group_layout.addWidget(pos_label)
            pos_layout = QHBoxLayout()
            for adduct in adducts_by_charge[charge]['positive']:
                cb = QCheckBox(adduct)
                cb.setChecked(True)
                self.adduct_checkboxes_by_charge[charge][adduct] = cb
                pos_layout.addWidget(cb)
            pos_layout.addStretch()
            group_layout.addLayout(pos_layout)

            # Negative adducts
            neg_label = QLabel("Negative Mode:")
            neg_label.setStyleSheet("font-style: italic; color: #cc6600;")
            group_layout.addWidget(neg_label)
            neg_layout = QHBoxLayout()
            for adduct in adducts_by_charge[charge]['negative']:
                cb = QCheckBox(adduct)
                cb.setChecked(True)
                self.adduct_checkboxes_by_charge[charge][adduct] = cb
                neg_layout.addWidget(cb)
            neg_layout.addStretch()
            group_layout.addLayout(neg_layout)

            group_box.setLayout(group_layout)
            self.adduct_group_boxes[charge] = group_box
            selection_layout.addWidget(group_box)

        selection_group.setLayout(selection_layout)

        # Initialize visibility based on default charge state selection
        self.update_adduct_visibility()

        # === ISOTOPE LABELING OPTIONS ===
        label_group = QGroupBox("Isotope Labeling Options")
        label_layout = QVBoxLayout()

        self.add_labels_checkbox = QCheckBox("Add isotope labels")
        label_layout.addWidget(self.add_labels_checkbox)

        isotope_layout = QHBoxLayout()
        isotope_layout.addWidget(QLabel("GSL Isotope:"))
        self.isotope_input = QLineEdit("M2DN15")
        isotope_layout.addWidget(self.isotope_input)
        label_layout.addLayout(isotope_layout)

        cer_layout = QHBoxLayout()
        cer_layout.addWidget(QLabel("Cer Isotope:"))
        self.cer_isotope_input = QLineEdit("M2DN15")
        cer_layout.addWidget(self.cer_isotope_input)
        label_layout.addLayout(cer_layout)

        doxcer_layout = QHBoxLayout()
        doxcer_layout.addWidget(QLabel("doxCer Isotope:"))
        self.doxcer_isotope_input = QLineEdit("M3D")
        doxcer_layout.addWidget(self.doxcer_isotope_input)
        label_layout.addLayout(doxcer_layout)

        lcb_layout = QHBoxLayout()
        lcb_layout.addWidget(QLabel("Label Keywords:"))
        self.lcb_input = QLineEdit("LCB,precursor,HG(-Hex")
        lcb_layout.addWidget(self.lcb_input)
        label_layout.addLayout(lcb_layout)

        self.blank_mz_checkbox = QCheckBox("Blank all m/z values")
        label_layout.addWidget(self.blank_mz_checkbox)

        label_group.setLayout(label_layout)

        # === SIDE-BY-SIDE LAYOUT ===
        horizontal_container = QHBoxLayout()
        horizontal_container.addWidget(selection_group)
        horizontal_container.addWidget(label_group)
        layout.addLayout(horizontal_container)

        # Buttons
        btn_layout = QHBoxLayout()
        self.config_btn = QPushButton("Configure Lipid Chains")
        self.run_btn = QPushButton("Generate Transitions")
        self.save_btn = QPushButton("Save CSV")
        self.clear_btn = QPushButton("Clear Results")
        self.reset_btn = QPushButton("Reset to Defaults")
        self.about_btn = QPushButton("About")
        self.exit_btn = QPushButton("Exit")
        btn_layout.addWidget(self.config_btn)
        btn_layout.addWidget(self.run_btn)
        btn_layout.addWidget(self.save_btn)
        btn_layout.addWidget(self.clear_btn)
        btn_layout.addWidget(self.reset_btn)
        btn_layout.addWidget(self.about_btn)
        btn_layout.addWidget(self.exit_btn)
        layout.addLayout(btn_layout)

        # Status label
        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        # Results table
        layout.addWidget(QLabel("Results:"))
        self.table = QTableWidget()
        layout.addWidget(self.table)

        scroll.setWidget(main_widget)
        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll)
        self.setLayout(main_layout)

        # Connect signals
        self.config_btn.clicked.connect(self.show_config_dialog)
        self.run_btn.clicked.connect(self.run_gen)
        self.save_btn.clicked.connect(self.save_csv)
        self.clear_btn.clicked.connect(self.clear_results)
        self.reset_btn.clicked.connect(self.reset_to_defaults)
        self.about_btn.clicked.connect(self.show_about)
        self.exit_btn.clicked.connect(self.exit_application)

        self.df = None
        self.update_info()
        self.load_ui_state()
        self.connect_auto_save_signals()


    def update_info(self):
        """Update info display when lipid class changes"""
        lipid_class = self.lipid_combo.currentText()
        structure = gslgen_v10.LipidDatabase.get_structure_description(lipid_class)
        mw_range = gslgen_v10.LipidDatabase.molecular_weight_range(lipid_class)
        sialic_count = gslgen_v10.LipidDatabase.get_sialic_acid_count(lipid_class)

        info = f"Structure: {structure}\n"
        info += f"MW Range: {mw_range} Da\n"
        if sialic_count > 0:
            info += f"Sialic Acids: {sialic_count}"

        self.info_text.setText(info)

    def set_recommended_charges(self):
        """Set recommended charge states for current lipid class"""
        lipid_class = self.lipid_combo.currentText()
        recommended = gslgen_v10.get_recommended_charges_for_lipid(lipid_class)

        for charge, checkbox in self.charge_checkboxes.items():
            checkbox.setChecked(charge in recommended)

    def show_config_dialog(self):
        """Show configuration dialog"""
        dialog = ConfigDialog(self)
        dialog.exec()

    def run_gen(self):
        """Generate transitions"""
        lipid_class = self.lipid_combo.currentText()

        # Get selections using new methods
        charge_states = self.get_selected_charge_states()
        selected_adducts = self.get_selected_adducts()

        # Validation
        if not charge_states:
            QMessageBox.warning(self, "No Charge States",
                            "Please select at least one charge state!")
            return

        if not selected_adducts:
            QMessageBox.warning(self, "No Adducts",
                            "Please select at least one adduct!")
            return

        try:
            logger.info(f"Starting generation: {lipid_class}, charges={charge_states}")
            # Generate transitions with configuration
            transitions = gslgen_v10.generate_transitions(
                lipid_class,
                charge_states,
                selected_adducts
            )

            logger.info(f"Generation complete: {len(transitions)} transitions")

            df = pd.DataFrame(transitions)

            # Add isotope labels if requested
            if self.add_labels_checkbox.isChecked():
                df = gslgen_v10.add_isotope_labels(
                    df,
                    isotope=self.isotope_input.text(),
                    doxcer_isotope=self.doxcer_isotope_input.text(),
                    cer_isotope=self.cer_isotope_input.text(),
                    lcb=self.lcb_input.text()
                )

            # Apply blanking
            if self.blank_mz_checkbox.isChecked():
                df = gslgen_v10.blank_mz_values(df)

            self.df = df
            self.populate_table(df)

            # Update status label
            label_info = " with isotope labels" if self.add_labels_checkbox.isChecked() else ""
            self._show_success(f"Generated {len(df)} transitions for {lipid_class}")



        except Exception as e:
            logger.error(f"Generation failed: {str(e)}", exc_info=True)
            self._show_error(f"Error: {str(e)}")

    def populate_table(self, df):
        """Populate table from dataframe"""
        self.table.setRowCount(len(df))
        self.table.setColumnCount(len(df.columns))
        self.table.setHorizontalHeaderLabels(df.columns)

        for i, row in df.iterrows():
            for j, col in enumerate(df.columns):
                self.table.setItem(i, j, QTableWidgetItem(str(row[col])))

        self.table.resizeColumnsToContents()

    def clear_results(self):
        """Clear results table"""
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.df = None
        self.status_label.setText("")

    def save_csv(self):
        """Save results to CSV"""
        if self.df is None:
            QMessageBox.warning(self, "No Data", "Generate transitions first!")
            return

        lipid_class = self.lipid_combo.currentText()
        default_filename = f"{lipid_class.lower()}_transitions.csv"

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save CSV",
            default_filename,
            "CSV Files (*.csv)"
        )

        if path:
            try:
                self.df.to_csv(path, index=False)
                self._show_success(f"Saved {len(self.df)} rows to {path}")
            except Exception as e:
                self._show_error(f"Save error: {str(e)}")

    def show_about(self):
        """Show about dialog with authorship info"""
        citation_link = "https://github.com/ahuelsmeier/gsl-transition-generator"

        about_text = """
        <h2>GSL & Ceramide Transition Generator v1.0</h2>
        <p><b>Author:</b> Andreas J. Hülsmeier</p>
        <p><b>Organization:</b> University of Zurich, University Hospital Zurich</p>
        <p><b>Email:</b> andreas.huelsmeier@uzh.ch</p>
        <p><b>Year:</b> 2025</p>
        <hr>
        <p>A tool for generating mass spectrometry transition lists
        for glycosphingolipids and ceramides.</p>
        <p><b>License:</b> MIT License</p>
        <p><b>Citation:</b> If you use this software in your research,
please visit the project repository for citation information:</p>
        <p><a href="{0}">{0}</a></p>
        """.format(citation_link)

        QMessageBox.about(self, "About GSL Transition Generator", about_text)

    def set_window_icon(self):
        """Set window icon based on platform"""
        import platform
        import os

        system = platform.system()

        if system == "Windows":
            icon_file = "app_icon.ico"
        elif system == "Darwin":
            icon_file = "app_icon.icns"
        else:
            icon_file = "icon_256.png"

        if os.path.exists(icon_file):
            self.setWindowIcon(QIcon(icon_file))
        else:
            print(f"Warning: Icon file '{icon_file}' not found")

    def exit_application(self):
        """Exit the application with confirmation"""
        if self.df is not None:
            reply = QMessageBox.question(
                self,
                'Confirm Exit',
                'Are you sure you want to exit? Any unsaved data will be lost.',
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                QApplication.quit()
        else:
            QApplication.quit()

    def _show_success(self, message):
        """Display a success status message"""
        self.status_label.setText(message)
        self.status_label.setStyleSheet(
            "color: green; "
            "font-weight: bold; "
            "background-color: #e8f5e9; "
            "padding: 5px; "
            "border-radius: 3px;"
        )

    def _show_error(self, message):
        """Display an error status message"""
        self.status_label.setText(f"❌ {message}")
        self.status_label.setStyleSheet(
            "color: red; "
            "font-weight: bold; "
            "background-color: #ffebee; "
            "padding: 5px; "
            "border-radius: 3px;"
        )

    def _show_warning(self, message):
        """Display a warning status message"""
        self.status_label.setText(f"⚠️ {message}")
        self.status_label.setStyleSheet(
            "color: #ff9800; "
            "font-weight: bold; "
            "background-color: #fff3e0; "
            "padding: 5px; "
            "border-radius: 3px;"
        )

    def _clear_status(self):
        """Clear the status message"""
        self.status_label.setText("")
        self.status_label.setStyleSheet("")

    def update_adduct_visibility(self):
        """Show/hide adduct groups based on selected charge states."""
        for charge, group_box in self.adduct_group_boxes.items():
            group_box.setVisible(self.charge_checkboxes[charge].isChecked())

    def get_selected_charge_states(self):
        """Return list of selected charge states."""
        return [charge for charge, cb in self.charge_checkboxes.items() if cb.isChecked()]

    def get_selected_adducts(self):
        """Return list of selected adducts from visible charge state groups."""
        selected = []
        for charge in self.get_selected_charge_states():
            for adduct, cb in self.adduct_checkboxes_by_charge[charge].items():
                if cb.isChecked():
                    selected.append(adduct)
        return selected or None

    def save_ui_state(self):
        """Save current UI selections to config"""
        if self.block_auto_save:  # ← ADD THIS CHECK
            return
        config = gslgen_v10.ConfigManager.load_config()
        config['charge_states'] = self.get_selected_charge_states()
        config['selected_adducts'] = self.get_selected_adducts()
        config['isotope_labeling'] = {
            'enabled': self.add_labels_checkbox.isChecked(),
            'gsl_isotope': self.isotope_input.text(),
            'cer_isotope': self.cer_isotope_input.text(),
            'doxcer_isotope': self.doxcer_isotope_input.text(),
            'label_keywords': self.lcb_input.text(),
            'blank_mz': self.blank_mz_checkbox.isChecked()
        }
        gslgen_v10.ConfigManager.save_config(config)

    def load_ui_state(self):
        """Load UI selections from config"""
        config = gslgen_v10.ConfigManager.load_config()
        saved_charges = config.get('charge_states', [1])
        for charge, cb in self.charge_checkboxes.items():
            cb.setChecked(charge in saved_charges)

        saved_adducts = config.get('selected_adducts', None)
        if saved_adducts:
            for charge in self.adduct_checkboxes_by_charge:
                for adduct, cb in self.adduct_checkboxes_by_charge[charge].items():
                    cb.setChecked(adduct in saved_adducts)

        isotope_config = config.get('isotope_labeling', {})
        self.add_labels_checkbox.setChecked(isotope_config.get('enabled', False))
        self.isotope_input.setText(isotope_config.get('gsl_isotope', 'M2DN15'))
        self.cer_isotope_input.setText(isotope_config.get('cer_isotope', 'M2DN15'))
        self.doxcer_isotope_input.setText(isotope_config.get('doxcer_isotope', 'M3D'))
        self.lcb_input.setText(isotope_config.get('label_keywords', 'LCB,precursor,HG(-Hex'))
        self.blank_mz_checkbox.setChecked(isotope_config.get('blank_mz', False))
        self.update_adduct_visibility()

    def connect_auto_save_signals(self):
        """Connect signals for auto-saving UI state"""
        for cb in self.charge_checkboxes.values():
            cb.stateChanged.connect(self.save_ui_state)
        for charge in self.adduct_checkboxes_by_charge:
            for cb in self.adduct_checkboxes_by_charge[charge].values():
                cb.stateChanged.connect(self.save_ui_state)
        self.add_labels_checkbox.stateChanged.connect(self.save_ui_state)
        self.isotope_input.textChanged.connect(self.save_ui_state)
        self.cer_isotope_input.textChanged.connect(self.save_ui_state)
        self.doxcer_isotope_input.textChanged.connect(self.save_ui_state)
        self.lcb_input.textChanged.connect(self.save_ui_state)
        self.blank_mz_checkbox.stateChanged.connect(self.save_ui_state)

    def reset_to_defaults(self):
        """Reset all settings to default values"""
        reply = QMessageBox.question(
            self,
            "Reset to Defaults",
            "This will reset all settings to their default values.\n\nAre you sure?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            # Block auto-save during reset
            self.block_auto_save = True

            # Get default config and save it
            default_config = gslgen_v10.ConfigManager.get_default_config()
            gslgen_v10.ConfigManager.save_config(default_config)

            # Reload UI from defaults
            self.load_ui_state()

            # Re-enable auto-save
            self.block_auto_save = False

            QMessageBox.information(
                self,
                "Reset Complete",
                "All settings have been reset to default values."
            )


if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = GSLGui()
    win.show()
    sys.exit(app.exec())

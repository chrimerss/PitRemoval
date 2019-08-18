# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PitRemovalDialog
                                 A QGIS plugin
 this plugin removes pits caused by artifacts
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2019-08-15
        git sha              : $Format:%H$
        copyright            : (C) 2019 by Zhi
        email                : chrimerss@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from qgis.PyQt import uic
from qgis.PyQt import QtWidgets
from qgis.PyQt.QtWidgets import QAction,QFileDialog, QMessageBox
from qgis.core import QgsProject,QgsRasterLayer

# This loads your .ui file so that PyQt can populate your plugin with the elements from Qt Designer
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'pit_removal_dialog_base.ui'))


class PitRemovalDialog(QtWidgets.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(PitRemovalDialog, self).__init__(parent)
        # Set up the user interface from Designer through FORM_CLASS.
        # After self.setupUi() you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.add_rasters()
        self.add_mode()

    def browseOutputFile(self):
        filename = QFileDialog.getSaveFileName(self, "Select output file ","", '*.tif')
        self.outputPlace.setText(filename[0])


    def browseInputFile(self):
        filename= QFileDialog.getOpenFileName(self, "Select input raster file ","", '*.tif')
        self.selectRaster.addItem(filename[0])
        self.selectRaster.setItemText(0,str(filename[0]))

    def add_rasters(self):
        root= QgsProject.instance().layerTreeRoot()
        layers_name= []
        for Layer in root.children():
            #check if raster
            if isinstance(Layer.layer(),QgsRasterLayer) and Layer.layer() not in root.children():
                self.selectRaster.addItem(Layer.name())

    def add_mode(self):
        mode= ['clip', 'minimum cost', 'balancing cost']
        for i in range(3):
            self.selectMode.addItem(mode[i])

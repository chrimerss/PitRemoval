# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PitRemoval
                                 A QGIS plugin
 this plugin removes pits caused by artifacts
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2019-08-15
        copyright            : (C) 2019 by Zhi
        email                : chrimerss@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load PitRemoval class from file PitRemoval.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .pit_removal import PitRemoval
    return PitRemoval(iface)

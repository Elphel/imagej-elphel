package com.elphel.imagej.common;

import java.util.Properties;

public class EProperties extends Properties{
	private static final long serialVersionUID = -425120416815883045L;
	public int getProperty(String key, int value){
		return Integer.parseInt(getProperty(key, ""+value));
	}
	public double getProperty(String key, double value){
		return Double.parseDouble(getProperty(key, ""+value));
	}
	public boolean getProperty(String key, boolean value){
		return Boolean.parseBoolean(getProperty(key, ""+value));
	}
}

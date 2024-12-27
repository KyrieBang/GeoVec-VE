package com.geosparkviz

import scala.math.pow
import scala.util.Random
import scala.io.Source
import scala.collection.mutable.ArrayBuffer
import org.apache.sedona.spark.SedonaContext
import org.apache.sedona.viz.core.Serde.SedonaVizKryoRegistrator
import org.apache.sedona.viz.sql.utils.SedonaVizRegistrator
import org.apache.sedona.common.enums.FileDataSplitter
import org.apache.sedona.core.formatMapper.shapefileParser.ShapefileReader
import org.apache.sedona.core.spatialRDD.{PointRDD, LineStringRDD, PolygonRDD, RectangleRDD, SpatialRDD}
import org.apache.sedona.viz.core.{ImageGenerator, ImageSerializableWrapper}
import org.apache.sedona.viz.extension.visualizationEffect.{ScatterPlot}
import org.apache.sedona.viz.utils.ImageType
import org.apache.spark.sql.SparkSession
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.Geometry
import java.awt.Color
import java.nio.file.{Files, Paths, StandardOpenOption}



object GeoVecViz {
  def readData(sedona: SparkSession, data_path: String, data_type: String): SpatialRDD[_] = {
    if (data_type == "point"){
      var data_RDD = ShapefileReader.readToPointRDD(sedona.sparkContext, data_path)
      return data_RDD
    }
    else if (data_type == "linestring"){
      var data_RDD = ShapefileReader.readToLineStringRDD(sedona.sparkContext, data_path)
      return data_RDD
    }
    else if (data_type == "polygon"){
      var data_RDD = ShapefileReader.readToPolygonRDD(sedona.sparkContext, data_path)
      return data_RDD
    }
    else{
      val emptyRDD = sedona.sparkContext.emptyRDD[Geometry]
      val spatialRDD = new SpatialRDD[Geometry]()
      spatialRDD.setRawSpatialRDD(emptyRDD)
      return spatialRDD
    }
  }

  def calTime(sedona: SparkSession, data_RDD: SpatialRDD[_]): Unit = {
    val s_time = System.currentTimeMillis()

    val viz_box = new Envelope(-180.0, 180.0, -90.0, 90.0)
    val imageGenerator = new ImageGenerator
    val temp_dir = Files.createTempDirectory("tmp-pic")
    val dir_path = temp_dir.toFile.getAbsolutePath

    for(level <-0 to 8){
      val slice_count = Math.pow(2, level)
      val total_width = 256 * Math.pow(2, level)
      val visual_operator = new ScatterPlot(total_width.toInt, total_width.toInt, viz_box, false, slice_count.toInt, slice_count.toInt, true, true)
      visual_operator.CustomizeColor(255, 255, 255, 255, Color.GREEN, true)
      visual_operator.Visualize(sedona.sparkContext, data_RDD)
      imageGenerator.SaveRasterImageAsLocalFile(visual_operator.distributedRasterImage, dir_path, ImageType.PNG)
    }

    val e_time = System.currentTimeMillis()
    val time_use: Double = (e_time - s_time) / 1000
    System.out.println(time_use)

    Files.walk(temp_dir)
      .sorted(java.util.Comparator.reverseOrder())
      .forEach(Files.delete)
  }

  def calSpeed(sedona: SparkSession, data_RDD: SpatialRDD[_]): Unit = {
    val viz_box = new Envelope(-180.0, 180.0, -90.0, 90.0)
    val imageGenerator = new ImageGenerator
    val dir_path = "../project/GeoSparkViz/tmp/"
    
    for(level <- 1 to 9 by 2){
      val slice_count = Math.pow(2, level)
      val total_width = 256 * Math.pow(2, level)
      
      val s_time = System.currentTimeMillis()

      val visual_operator = new ScatterPlot(total_width.toInt, total_width.toInt, viz_box, false, slice_count.toInt, slice_count.toInt, true, true)
      visual_operator.CustomizeColor(255, 255, 255, 255, Color.GREEN, true)
      
      visual_operator.Visualize(sedona.sparkContext, data_RDD)
      imageGenerator.SaveRasterImageAsLocalFile(visual_operator.distributedRasterImage, dir_path, ImageType.PNG)
      
      val e_time = System.currentTimeMillis()
      val time_use: Double = (e_time - s_time) / 1000
      
      val pic_count = Files.list(Paths.get(dir_path)).count()
      System.out.println(pic_count / time_use)

      Files.walk(Paths.get(dir_path))
      .filter(Files.isRegularFile(_))
      .forEach(Files.delete(_))
    }
  }


  def main(args: Array[String]): Unit = {
    // define path
    val data_path = args(0)
    val data_type = args(1)
    val plot_type = args(2)
    
    // Initialize Spark context
    val config = SedonaContext.builder().appName("GeoSparkViz")
      .master("local[*]")
      .config("spark.kryo.registrator", classOf[SedonaVizKryoRegistrator].getName)
      .getOrCreate()
    val sedona = SedonaContext.create(config)
    SedonaVizRegistrator.registerAll(sedona)

    var data_RDD = readData(sedona, data_path, data_type)

    if (plot_type == "fig15"){
      calTime(sedona, data_RDD)
    }
    else if (plot_type == "fig16"){
      calSpeed(sedona, data_RDD)
    }
  }
}


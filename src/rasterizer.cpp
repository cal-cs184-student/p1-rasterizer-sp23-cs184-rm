#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
      
    int s_width = width * (sqrt(sample_rate));
    sample_buffer[y * s_width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // determines if the point is inside the triange
  int pointinside(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {
      float L0 = -1 * (x - x0) * (y1 - y0) + (y - y0) * (x1 - x0);
      float L1 = -1 * (x - x1) * (y2 - y1) + (y - y1) * (x2 - x1);
      float L2 = -1 * (x - x2) * (y0 - y2) + (y - y2) * (x0 - x2);
      
      if ((L0 > 0) && (L1 > 0) && (L2 > 0)) {
          return 1;
      } else if ((L0 < 0) && (L1 < 0) && (L2 < 0)) {
          return 1;
      }
      return 0;
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    int xmax = (int) max(max(x0, x1), x2);
    int xmin = (int) min(min(x0, x1), x2);
    int ymax = (int) max(max(y0, y1), y2);
    int ymin = (int) min(min(y0, y1), y2);
//    for (int x = xmin; x < xmax; x++) {
//        for (int y = ymin; y < ymax; y++) {
//            if (pointinside(x + .5, y + .5, x0, y0, x1, y1, x2, y2) == 1) {
//                fill_pixel(x, y, color);
//            }
//        }
//    }
    // TODO: Task 2: Update to implement super-sampled rasterization
      
      int srs = (int) sqrt(sample_rate);
      float xpix;
      float ypix;
      float s_width = width * srs;
      for (int x = xmin; x < xmax; x++) {
          for (int y = ymin; y < ymax; y++) {
              for (int i = 0; i < srs; i++) {
                  for (int j = 0; j < srs; j++) {
                      xpix = x + ((i + 1) / srs + i / srs) / 2;
                      ypix = y + ((j + 1) / srs + j / srs) / 2;
                      if (pointinside(xpix, ypix, x0, y0, x1, y1, x2, y2) == 1) {
                          fill_pixel(srs * x + i, srs * y + j, color);
                      }
                  }
              }
          }
      }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      int xmax = (int) max(max(x0, x1), x2);
      int xmin = (int) min(min(x0, x1), x2);
      int ymax = (int) max(max(y0, y1), y2);
      int ymin = (int) min(min(y0, y1), y2);
  //    int image[xmax][ymax];
      for (int x = xmin; x < xmax; x++) {
          for (int y = ymin; y < ymax; y++) {
  //            fill_pixel(x, y, color);
              if (pointinside(x + .5, y + .5, x0, y0, x1, y1, x2, y2) == 1) {
                  float alpha = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                  float beta = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                  float gamma = 1 - alpha - beta;
                  float red = alpha * c0.r + beta * c1.r + gamma * c2.r;
                  float green = alpha * c0.g + beta * c1.g + gamma * c2.g;
                  float blue = alpha * c0.b + beta * c1.b + gamma * c2.b;
                  fill_pixel(x, y, Color(red, green, blue));
              }
          }
      }


  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
      int xmax = (int) max(max(x0, x1), x2);
      int xmin = (int) min(min(x0, x1), x2);
      int ymax = (int) max(max(y0, y1), y2);
      int ymin = (int) min(min(y0, y1), y2);
      for (int x = xmin; x < xmax; x++) {
          for (int y = ymin; y < ymax; y++) {
              SampleParams sp;
              sp.lsm = lsm;
              sp.psm = psm;
              sp.p_uv = Vector2D(u0 / width, v0 / height);
              Color col = tex.sample(sp);
          }
      }
      




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = sqrt(rate) * sqrt(rate);


    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
      
    float srs = sqrt(sample_rate);
    float s_width = width * srs;
    float red = 0;
    float blue = 0;
    float green = 0;
;
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
//        Color col = sample_buffer[y * width + x];
          for (int i = 0; i < srs; i++) {
              for (int j = 0; j < srs; j++) {
                  Color c = sample_buffer[(srs * y + j) * s_width + (srs * x + i)];
                  red += c.r;
                  blue += c.b;
                  green += c.g;
              }
          }
          red = red / sample_rate;
          blue = blue / sample_rate;
          green = green / sample_rate;
          Color col = Color(red, blue, green);

          for (int k = 0; k < 3; ++k) {
            this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
          }
          red = 0;
          blue = 0;
          green = 0;
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL

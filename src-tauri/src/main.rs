// target_os: linux, macos, windows
#![cfg_attr(
  all(not(debug_assertions), target_os = "macos"),
)]

use std::sync::Arc;

use tauri::async_runtime::RwLock;

use std::borrow::BorrowMut;
mod diffusion;
mod menu;
mod simulation;
mod staggered_grid;
mod utils;

#[derive(Debug, Clone)]
struct AppState{
	pub info: utils::FluidInformation,
  pub grid: Option<staggered_grid::StaggeredGrid>,
}

impl AppState {
  fn new() -> Self {
    Self {
      info: utils::FluidInformation::new(),
      grid: None,
    }
  }
  pub fn update_grid(&mut self, width: usize, height: usize) {
    self.grid = Some(staggered_grid::StaggeredGrid::new(height, width));
  }
}
#[derive(Clone)]
pub struct ConcurrentAppState {
  inner: Arc<RwLock<AppState>>,
}

impl ConcurrentAppState {
  pub fn new() -> Self {
    Self {
      inner: Arc::new(RwLock::new(AppState::new())),
    }
  }
  pub async fn update_grid_size(&self, width: usize, height: usize) {
    let mut lock = self.inner.write().await;
    lock.grid = Some(staggered_grid::StaggeredGrid::new(width, height));
    lock.info.width = width;
    lock.info.height = height;
  }
}
fn main() {
  println!("Hello, world!");
  tauri::Builder::default()
    .manage(ConcurrentAppState::new())
    .invoke_handler(tauri::generate_handler![
      initialize_state,
      get_grid_colors,
      take_step,
    ])
    .menu(menu::mainmenu())
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}

#[tauri::command]
async fn initialize_state(width: usize, height: usize, state: tauri::State<'_, ConcurrentAppState>) -> Result<(), String> {
  println!("initialize_state");
  if state.inner.read().await.grid.is_none() {
    return Err("Grid already initialized".to_string());
  }
  state.update_grid_size(width, height).await;
  println!("initialize_state done");
  Ok(())
}
#[tauri::command]
async fn get_grid_colors(state: tauri::State<'_, ConcurrentAppState>) -> Result<Vec<u8>, String> {
  //TODO: Refactor this to be a state variable
  let positive_color = utils::Color {r: 1., g: 0., b: 0.}; // red
  let negative_color = utils::Color {r: 0., g: 178./255., b: 1.}; // blue
  let lock = state.inner.read().await;
  if lock.grid.is_none() {
    return Err("Grid not initialized".to_string());
  }
  let grid = lock.grid.as_ref().unwrap();
  let shape = grid.contents.shape();
  let height = shape[0];
  let width = shape[1];
  let mut colors: Vec<u8> = Vec::with_capacity(height * width * 4);
  println!("adding colors with shape {:?}", grid.contents.shape());
  for row in 0..height {
    for col in 0..width {
      let val = grid.contents[(row, col)].clamp(-1., 1.);
      let color = utils::lerp_color(val, negative_color, positive_color);
      let color_index = row * width + col;
      colors[color_index] = color.r as u8 * 255;
      colors[color_index + 1] = color.g as u8 * 255;
      colors[color_index + 2] = color.b as u8 * 255;
      colors[color_index + 3] = 255;
    }
  }
  println!("colors added");
  Ok(colors)
}

#[tauri::command]
async fn take_step(state: tauri::State<'_, ConcurrentAppState>) -> Result<(), String> {
  let lock = state.inner.read().await;
  if lock.grid.is_none() {
    return Err("Grid not initialized".to_string());
  }
  // TODO: Refactor this to not use clone
  let mut grid = lock.grid.as_ref().unwrap().clone();
  let info = lock.info.clone();
  simulation::simulate_step(&mut grid, 1.0, &info);
  Ok(())
}
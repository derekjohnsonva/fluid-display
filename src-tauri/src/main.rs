// target_os: linux, macos, windows
#![cfg_attr(all(not(debug_assertions), target_os = "macos"),)]
use std::sync::Arc;

use tauri::async_runtime::RwLock;

mod diffusion;
mod menu;
mod simulation;
mod staggered_grid;
mod utils;

#[derive(Debug, Clone)]
struct AppState {
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

  pub async fn simulate_step(&self) -> Result<(), String> {
    //TODO: This should be done without cloning
    let mut grid = self.inner.read().await.grid.as_ref().unwrap().clone();
    simulation::simulate_step(&mut grid, 1.0, &self.inner.read().await.info);
    self.inner.write().await.grid = Some(grid);
    Ok(())
  }
  pub async fn update_heat_transfer(&self, new_ht_rate: f32) {
    self.inner.write().await.info.heat_transfer_rate = new_ht_rate;
  }
  pub async fn update_bouyancy(&self, new_bouyancy: f32) {
    self.inner.write().await.info.bouyancy = new_bouyancy;
  }
  pub async fn update_diffusion(&self, new_diffusion: f32) {
    if new_diffusion > 0. {
      self.inner.write().await.info.diffuse = Some(new_diffusion);
    } else {
      self.inner.write().await.info.diffuse = None;
    }
  }
  pub async fn update_viscosity(&self, new_viscosity: f32) {
    if new_viscosity > 0. {
      self.inner.write().await.info.viscosity = Some(new_viscosity);
    } else {
      self.inner.write().await.info.viscosity = None;
    }
  }
}
fn main() {
  tauri::Builder::default()
    .manage(ConcurrentAppState::new())
    .invoke_handler(tauri::generate_handler![
      initialize_state,
      get_grid_colors,
      take_step,
      update_heat_transfer,
      update_bouyancy,
      update_diffusion,
      update_viscosity,
    ])
    .menu(menu::mainmenu())
    .run(tauri::generate_context!())
    .expect("error while running tauri application");
}

#[tauri::command]
async fn initialize_state(
  width: usize,
  height: usize,
  state: tauri::State<'_, ConcurrentAppState>,
) -> Result<(), String> {
  if state.inner.read().await.grid.is_some() {
    return Err("Grid already initialized".to_string());
  }
  state.update_grid_size(width, height).await;
  Ok(())
}
#[tauri::command]
async fn get_grid_colors(state: tauri::State<'_, ConcurrentAppState>) -> Result<Vec<u8>, String> {
  //TODO: Refactor this to be a state variable
  // println!("starting grid colors");
  let positive_color = utils::Color {
    r: 1.,
    g: 0.,
    b: 0.,
  }; // red
  let negative_color = utils::Color {
    r: 0.,
    g: 178. / 255.,
    b: 1.,
  }; // blue
  let lock = state.inner.read().await;
  if lock.grid.is_none() {
    return Err("Grid not initialized while trying to get grid colors".to_string());
  }
  let grid = lock.grid.as_ref().unwrap();
  let shape = grid.contents.shape();
  let height = shape[0];
  let width = shape[1];
  let mut colors: Vec<u8> = Vec::with_capacity(height * width * 4);
  // println!("adding colors with shape {:?}", grid.contents.shape());
  for row in 0..height {
    for col in 0..width {
      let val = grid.contents[(row, col)].clamp(-1., 1.);
      let color = utils::lerp_color(val, negative_color, positive_color);
      colors.push((color.r * 255.) as u8);
      colors.push((color.g * 255.) as u8);
      colors.push((color.b * 255.) as u8);
      colors.push(255);
    }
  }
  println!("got grid colors");
  Ok(colors)
}

#[tauri::command]
async fn take_step(state: tauri::State<'_, ConcurrentAppState>) -> Result<(), String> {
  state.simulate_step().await;
  println!("took step");
  Ok(())
}

#[tauri::command]
async fn update_heat_transfer(
  state: tauri::State<'_, ConcurrentAppState>,
  new_val: f32,
) -> Result<(), String> {
  state.update_heat_transfer(new_val).await;
  println!("updated heat transfer");
  Ok(())
}
#[tauri::command]
async fn update_bouyancy(
  state: tauri::State<'_, ConcurrentAppState>,
  new_val: f32,
) -> Result<(), String> {
  state.update_bouyancy(new_val).await;
  println!("updated heat transfer");
  Ok(())
}
#[tauri::command]
async fn update_diffusion(
  state: tauri::State<'_, ConcurrentAppState>,
  new_val: f32,
) -> Result<(), String> {
  state.update_diffusion(new_val).await;
  println!("updated heat transfer");
  Ok(())
}
#[tauri::command]
async fn update_viscosity(
  state: tauri::State<'_, ConcurrentAppState>,
  new_val: f32,
) -> Result<(), String> {
  state.update_viscosity(new_val).await;
  println!("updated heat transfer");
  Ok(())
}
